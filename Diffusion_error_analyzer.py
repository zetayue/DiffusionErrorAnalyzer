# coding: utf-8

from __future__ import division, unicode_literals

import numpy as np
import math

from pymatgen.core import structure, get_el_sp
import pymatgen.core.physical_constants as phyc
from monty.json import MSONable
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.util.coord_utils import pbc_diff


class DiffusionErrorAnalyzer(MSONable):

    def __init__(self, material, displace, specie, temp,
                 steptime, skipstep, smoothed = "constant", min_obs = 30,
                 avg_nsteps = 1000, mint = 1000, maxt = 0.1):

        self.material = material
        self.disp = displace
        self.specie = specie
        self.temp = temp
        self.steptime = steptime
        self.skipstep = skipstep
        self.min_obs = min_obs
        self.smoothed = smoothed
        self.avg_nsteps = avg_nsteps
        self.mint = mint
        self.maxt = maxt

        indices = []
        framework_indices = []
        self.ion_num = 0
		
        for i, site in enumerate(material):
            if site.specie.symbol == specie:
                indices.append(i)
                self.ion_num = self.ion_num + 1
            else:
                framework_indices.append(i)
        print 'Number of diffusion ions:',self.ion_num
		
        if self.disp.shape[1] < 2:
            self.diffus = 0.
            self.conduct = 0.
            self.diffus_comp = np.array([0., 0., 0.])
            self.conduct_comp = np.array([0., 0., 0.])
            self.max_disp = 0
        else:
            framework_disp = self.disp[framework_indices]
            drift = np.average(framework_disp, axis = 0)[None, :, :]
            dc = self.disp - drift
            df = material.lattice.get_fractional_coords(dc)

            nions, nsteps, dim = dc.shape

            if not smoothed:
                timesteps = np.arange(0, nsteps)
            elif smoothed == "constant":
                if nsteps <= avg_nsteps:
                    raise ValueError('Not enough data to calculate diffus')
                timesteps = np.arange(0, nsteps - avg_nsteps)
            else:
                min_dt = int(mint / (self.skipstep * self.steptime))
                max_dt = int(1 + (nsteps - mint / (self.skipstep * self.steptime))*self.maxt\
                          + mint / (self.skipstep * self.steptime))
                print 'Sequence number of first MD data:', min_dt
                print 'Sequence number of last MD data:', max_dt
                if min_dt >= max_dt:
                    raise ValueError('Not enough data to calculate diffus')
                timesteps = np.arange(min_dt, max_dt, 1)

            dt = timesteps * self.steptime * self.skipstep

            msd = np.zeros_like(dt, dtype=np.double)
            sq_disp_ions = np.zeros((len(dc), len(dt)), dtype=np.double)
            msd_comp = np.zeros(dt.shape + (3,))

            lengths = np.array(self.material.lattice.abc)[None, None, :]

            for i, n in enumerate(timesteps):
                if not smoothed:
                    dx = dc[:, i:i + 1, :]
                    dcomp = df[:, i:i + 1, :] * lengths
                elif smoothed == "constant":
                    dx = dc[:, i:i + avg_nsteps, :] - dc[:, 0:avg_nsteps, :]
                    dcomp = (df[:, i:i + avg_nsteps, :]
                                   - df[:, 0:avg_nsteps, :]) * lengths
                else:
                    dx = dc[:, n:, :] - dc[:, :-n, :]
                    dcomp = (df[:, n:, :] - df[:, :-n, :]) * lengths
                sq_disp = dx ** 2
                sq_disp_ions[:, i] = np.average(np.sum(sq_disp, axis=2), axis=1)
                msd[i] = np.average(sq_disp_ions[:, i][indices])

                msd_comp[i] = np.average(dcomp[indices] ** 2, axis=(0, 1))

            def weighted_square(a, b):
                return np.linalg.lstsq(a, b)

            m_comp = np.zeros(3)
            m_comp_res = np.zeros(3)
            a = np.ones((len(dt), 2))
            a[:, 0] = dt
            for i in range(3):
                (m, c), res, rank, s = weighted_square(a, msd_comp[:, i])
                m_comp[i] = max(m, 1e-15)
                m_comp_res[i] = res[0]

            (m, c), res, rank, s = weighted_square(a, msd)

            m = max(m, 1e-15)

            conv_factor = get_conversion_factor(self.material, self.specie, self.temp)
            self.diffus = m / 60
            n = len(dt)
            denom = (n * np.sum(dt ** 2) - np.sum(dt) ** 2) * (n - 2)

            self.diffus_std_dev = np.sqrt(n * res[0] / denom) / 60
            self.conduct = self.diffus * conv_factor
            self.conduct_std_dev = self.diffus_std_dev * conv_factor

            self.diffus_comp = m_comp / 20
            self.diffus_comp_std_dev = np.sqrt(n * m_comp_res / denom) / 20
            self.conduct_comp = self.diffus_comp * conv_factor
            self.conduct_comp_std_dev = self.diffus_comp_std_dev * conv_factor

            self.drift = drift
            self.corrected_displace = dc
            self.max_ion_displace = np.max(np.sum(dc ** 2, axis=-1) ** 0.5, axis=1)
            self.max_disp = np.max(self.max_ion_displace[framework_indices])

            self.msd = msd
            self.msd_true = (msd[-1] - msd[1]) * self.ion_num / self.maxt  
            self.jump = (msd[-1] - msd[1]) * self.ion_num / self.maxt / 4  
            self.sq_disp_ions = sq_disp_ions
            self.msd_comp = msd_comp
            self.dt = dt
            self.indices = indices
            self.framework_indices = framework_indices

            cutoff_time = 0
            for i in msd:
                cutoff_time=cutoff_time + 1
                if i > 10:
                    break
            print cutoff_time


    def get_results(self, include_msd_t=False):

        A = 2.11131 * 2 * math.pow((self.maxt * 100 - 1.80496), 0.04733)
        B = 0.03026 + 0.000566396 * self.maxt * 100
        print A
        print B
        self.msd_true = max(self.msd_true, 1e-15)

        if 1 - (A/math.sqrt(self.msd_true) + B) > 0:
            lowerlimit = math.log(self.diffus * (1 - \
			                   (A / math.sqrt(self.msd_true) + B))) / math.log(10)
        else:
            lowerlimit = 1e+10
        d = {
            "D": self.diffus,
            "LogD": math.log(self.diffus) / math.log(10),
            "D_sigma": self.diffus * (A / math.sqrt(self.msd_true) + B),
            "LogD_upper_limit": math.log(self.diffus * (1 + \
			                   (A / math.sqrt(self.msd_true) + B))) / math.log(10),
            "LogD_lower_limit": lowerlimit,
            "Jump_per_atom": self.jump,
            "S": self.conduct,
            "specie": str(self.specie),
            "skipstep": self.skipstep,
            "steptime": self.steptime,
            "temp": self.temp,
            "max_disp": self.max_disp
        }
        if include_msd_t:
            d["msd"] = self.msd.tolist()
            d["msd_true"] = self.msd_true.tolist()
            d["msd_comp"] = self.msd_comp.tolist()
            d["dt"] = self.dt.tolist()
        return d

    def get_plot(self, plt=None, mode="specie", element="Na"):

        from pymatgen.util.plotting_utils import get_publication_quality_plot
        plt = get_publication_quality_plot(12, 8, plt=plt)

        if mode == "species":
            for sp in sorted(self.material.composition.keys()):
                indices = [i for i, site in enumerate(self.material) if site.specie == sp]
                sd = np.average(self.sq_disp_ions[indices, :], axis=0)
                plt.plot(self.dt, sd, label=sp.__str__())
            plt.legend(loc=2, prop={"size": 20})
        elif mode == "ions":
            for i, site in enumerate(self.material):
                sd = self.sq_disp_ions[i, :]
                if site.specie.__str__() == element:
                    plt.plot(self.dt, sd)
            plt.legend(loc=2, prop={"size": 20})
        else: 
            plt.plot(self.dt, self.msd, 'k')
            plt.plot(self.dt, self.msd_comp[:, 0], 'r')
            plt.plot(self.dt, self.msd_comp[:, 1], 'g')
            plt.plot(self.dt, self.msd_comp[:, 2], 'b')
            plt.legend(["Overall", "a", "b", "c"], loc=2, prop={"size": 20})
        plt.xlabel("Timestep (fs)")
        plt.ylabel("MSD ($\AA^2$)")
        plt.tight_layout()
        return plt

    def plot(self, mode="default"):

        self.get_plot(mode=mode).show()

    @classmethod
    def from_structure(cls, structure, specie, temp,
                        steptime, skipstep, smoothed="constant", min_obs=30,
                        avg_nsteps=1000, initial_disp=None,
                        initial_material=None,mint=0.01, maxt=0.1):
        p = []
        for i, s in enumerate(structure):
            if i == 0:
                material = s
            p.append(np.array(s.frac_coords)[:, None])

        if initial_material is not None:
            p.insert(0, np.array(initial_material.frac_coords)[:, None])
        else:
            p.insert(0, p[0])
        p = np.concatenate(p, axis = 1)
        dp = p[:, 1:] - p[:, :-1]
        dp = dp - np.round(dp)
        f_disp = np.cumsum(dp, axis = 1)
        if initial_disp is not None:
            f_disp += structure.lattice.get_fractional_coords(initial_disp)[:,
                                                              None, :]
        disp = structure.lattice.get_cartesian_coords(f_disp)

        return cls(material, disp, specie, temp,
                   steptime, skipstep=skipstep, smoothed=smoothed,
                   min_obs=min_obs, avg_nsteps=avg_nsteps,mint=mint, maxt=maxt)

    @classmethod
    def from_vaspruns(cls, vaspruns, specie, smoothed = "constant", min_obs = 30,
                      avg_nsteps = 1000, initial_disp = None,
                      initial_material = None, mint = 0.01, maxt = 0.1):

        def get_structure(vaspruns):
            for i, vr in enumerate(vaspruns):
                if i == 0:
                    skipstep = vr.ionic_skipstep or 1
                    final_material = vr.initial_material
                    temp = vr.parameters['TEEND']
                    steptime = vr.parameters['POTIM']
                    yield skipstep, temp, steptime

                fdist = pbc_diff(vr.initial_material.frac_coords, final_material.frac_coords)
                if np.any(fdist > 0.001):
                    raise ValueError('initial and final structure do not ''match.')
                final_material = vr.final_material

                assert (vr.ionic_skipstep or 1) == skipstep
                for s in vr.ionic_steps:
                    yield s['material']

        s = get_structure(vaspruns)
        skipstep, temp, steptime = next(s)

        return cls.from_structure(structure=s, specie=specie,
            temp=temp, steptime=steptime, skipstep=skipstep,
            smoothed=smoothed, min_obs=min_obs, avg_nsteps=avg_nsteps,
            initial_disp=initial_disp, initial_material=initial_material,mint=mint, maxt=maxt)

    @classmethod
    def from_files(cls, filepaths, specie, skipstep = 40, smoothed = "constant",
                   min_obs = 30, avg_nsteps = 1000, ncores = None, initial_disp = None,
                   initial_material = None, mint = 0.01, maxt = 0.1):
        
        if ncores is not None and len(filepaths) > 1:
            import multiprocessing
            p = multiprocessing.Pool(ncores)
            vaspruns = p.imap(_get_vasprun,
                             [(fp, skipstep) for fp in filepaths])
            analyzer = cls.from_vaspruns(vaspruns, min_obs=min_obs,
                smoothed=smoothed, specie=specie, initial_disp=initial_disp,
                initial_material=initial_material, avg_nsteps=avg_nsteps)
            p.close()
            p.join()
            return analyzer
        else:
            def vr(filepaths):
                offset = 0
                for p in filepaths:
                    v = Vasprun(p, ionic_step_offset=offset, ionic_skipstep=skipstep)
                    yield v

                    offset = (-(v.nionic_steps - offset)) % skipstep
            return cls.from_vaspruns(vr(filepaths), min_obs=min_obs,
                smoothed=smoothed, specie=specie, initial_disp=initial_disp,
                initial_material=initial_material, avg_nsteps=avg_nsteps, mint=mint, maxt=maxt)


def _get_vasprun(args):
    
    return Vasprun(args[0], ionic_skipstep=args[1], parse_dos=False, parse_eigen=False)