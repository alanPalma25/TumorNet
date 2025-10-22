import numpy as np
import random, math, os, tempfile, time as pytime
import matplotlib.pyplot as plt
from matplotlib import colors
import imageio
import argparse
import configparser


# ======================================================
# === Parameter and Configuration System ===============
# ======================================================

class Parameters:
    """Holds all constants for the simulation, can load from CLI or .ini."""
    def __init__(self):
        # --- Default parameters ---
        self.Lx, self.Ly = 60, 60
        self.seed = 12

        # Cell types
        self.EMPTY, self.CSC, self.PC = 0, 1, 2

        # Rates
        self.rate_div_CSC = 0.12
        self.rate_div_PC = 0.20
        self.rate_death_CSC = 0.002
        self.rate_death_PC = 0.01
        self.rate_move = 0.04

        # CSC division outcomes
        self.p_sym_self = 0.12
        self.p_sym_diff = 0.04
        self.pc_max_divisions = 3

        # Nutrient PDE parameters
        self.use_nutrient = True
        self.D = 1.0
        self.decay = 0.01
        self.uptake = 0.02
        self.nutrient_source = 1.0
        self.nutrient_dt = 0.4

        # Therapy parameters
        self.use_therapy = True
        self.therapy_start = 10.0
        self.therapy_end = 18.0
        self.therapy_fold_increase_PC_death = 20.0
        self.therapy_fold_increase_CSC_death = 5.0

        # Simulation controls
        self.t_max = 2000.
        self.frames_count = 100
        self.snapshot_times = np.linspace(0, self.t_max, self.frames_count)

        # Output
        self.out_dir = 'output'
        self.out_gif = os.path.join(self.out_dir, 'tumor_dynamics.gif')
        self.out_csv = os.path.join(self.out_dir, 'simulation_data.csv')

    def apply_seed(self):
        np.random.seed(self.seed)
        random.seed(self.seed)

    # ---------- CONFIGURATION LOADING ---------------
    def load_from_ini(self, filename):
        """Load parameters from an .ini configuration file."""
        config = configparser.ConfigParser()
        config.read(filename)
        if 'SIMULATION' in config:
            for key, val in config['SIMULATION'].items():
                if hasattr(self, key):
                    casted_val = self._cast_value(val)
                    setattr(self, key, casted_val)

    def load_from_args(self, args):
        """Override parameters with command-line arguments (if provided)."""
        for key, val in vars(args).items():
            if val is not None and hasattr(self, key):
                setattr(self, key, val)

    def _cast_value(self, val):
        """Infer type of value from string."""
        for cast in (int, float):
            try: return cast(val)
            except ValueError: pass
        if val.lower() in ('true', 'false'):
            return val.lower() == 'true'
        return val


# ======================================================
# === Simulation Components ============================
# ======================================================

class Lattice:
    """Represents spatial lattice for cells and nutrients."""
    def __init__(self, params: Parameters):
        self.p = params
        self.cell_type = np.zeros((self.p.Lx, self.p.Ly), dtype=np.int8)
        self.pc_div_left = np.zeros((self.p.Lx, self.p.Ly), dtype=np.int8)
        self.nutrient = np.ones((self.p.Lx, self.p.Ly)) * self.p.nutrient_source if self.p.use_nutrient else None
        self.nbrs = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]

        cx, cy = self.p.Lx // 2, self.p.Ly // 3
        self.cell_type[cx-5, cy] = self.p.PC
        self.cell_type[cx+5, cy] = self.p.CSC
        self.pc_div_left[self.cell_type == self.p.PC] = self.p.pc_max_divisions

    def in_bounds(self, x, y): return 0 <= x < self.p.Lx and 0 <= y < self.p.Ly

    def get_empty_neighbors(self, x, y):
        coords = []
        for dx, dy in self.nbrs:
            nx, ny = x + dx, y + dy
            if self.in_bounds(nx, ny) and self.cell_type[nx, ny] == self.p.EMPTY:
                coords.append((nx, ny))
        return coords

    def nutrient_step(self, dt):
        if not self.p.use_nutrient: return
        new = self.nutrient.copy()
        lap = np.zeros_like(self.nutrient)
        lap[1:-1, 1:-1] = (
            self.nutrient[2:, 1:-1] + self.nutrient[:-2, 1:-1] +
            self.nutrient[1:-1, 2:] + self.nutrient[1:-1, :-2] -
            4 * self.nutrient[1:-1, 1:-1]
        )
        uptake_field = np.zeros_like(self.nutrient)
        uptake_field[self.cell_type != self.p.EMPTY] = self.p.uptake
        new += dt * (self.p.D * lap - self.p.decay * self.nutrient - uptake_field)
        new[0, :], new[-1, :], new[:, 0], new[:, -1] = [self.p.nutrient_source]*4
        new[new < 0] = 0.0
        self.nutrient = new


class EventSystem:
    """Handles Gillespie event calculations."""
    def __init__(self, lattice: Lattice):
        self.lattice = lattice
        self.p = lattice.p

    def compute_propensities(self, current_time):
        events, rates = [], []
        if self.p.use_therapy and self.p.therapy_start <= current_time <= self.p.therapy_end:
            rd_pc = self.p.rate_death_PC * self.p.therapy_fold_increase_PC_death
            rd_csc = self.p.rate_death_CSC * self.p.therapy_fold_increase_CSC_death
        else:
            rd_pc, rd_csc = self.p.rate_death_PC, self.p.rate_death_CSC

        xs, ys = np.nonzero(self.lattice.cell_type != self.p.EMPTY)
        for x, y in zip(xs, ys):
            c = self.lattice.cell_type[x, y]
            nut_factor = 1.0
            if self.p.use_nutrient:
                nloc = self.lattice.nutrient[x, y]
                nut_factor = nloc / (0.5 + nloc)

            if c == self.p.CSC:
                rdiv, rdeath, rmove = self.p.rate_div_CSC * nut_factor, rd_csc, self.p.rate_move
            elif c == self.p.PC:
                rdiv = self.p.rate_div_PC * nut_factor if self.lattice.pc_div_left[x, y] > 0 else 0.0
                rdeath, rmove = rd_pc, self.p.rate_move * 0.5
            else:
                continue

            if self.lattice.get_empty_neighbors(x, y):
                if rdiv > 0: events.append(('div', x, y)); rates.append(rdiv)
                events.append(('move', x, y)); rates.append(rmove)
            events.append(('death', x, y)); rates.append(rdeath)

        if not rates:
            return [], [], 0.0
        rates = np.array(rates, dtype=float)
        return events, rates, rates.sum()


class TumorSimulation:
    """Main simulation controller."""
    def __init__(self, params: Parameters):
        self.p = params
        self.p.apply_seed()
        self.lattice = Lattice(params)
        self.event_system = EventSystem(self.lattice)
        self.t, self.event_count = 0.0, 0
        self.max_events = 2_000_000
        self.tempdir = tempfile.mkdtemp()
        self.frames, self.frame_idx = [], 0
        self.next_snapshot_idx = 0
        self.data_records = []

    def record_state(self):
        csc = np.sum(self.lattice.cell_type == self.p.CSC)
        pc = np.sum(self.lattice.cell_type == self.p.PC)
        mean_n = np.mean(self.lattice.nutrient) if self.p.use_nutrient else 0
        self.data_records.append((self.t, csc, pc, mean_n))

    def capture_frame(self):
        cmap = colors.ListedColormap(['white', 'red', 'blue'])
        bounds = [0, 0.5, 1.5, 2.5]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        fig = plt.figure(figsize=(4, 4), dpi=100)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.imshow(self.lattice.cell_type.T, origin='lower', cmap=cmap, norm=norm)
        ax.text(2, 2, f"t={self.t:.2f}", color='black', bbox=dict(facecolor='white', alpha=0.6), fontsize=8)
        fname = os.path.join(self.tempdir, f"frame_{self.frame_idx:04d}.png")
        fig.savefig(fname, dpi=100)
        plt.close(fig)
        self.frames.append(fname)
        self.frame_idx += 1

    def run(self):
        start_wall = pytime.time()
        self.capture_frame()
        self.record_state()

        while self.t < self.p.t_max and self.event_count < self.max_events:
            events, rates, total_rate = self.event_system.compute_propensities(self.t)
            if total_rate <= 0: break
            dt = -math.log(np.random.rand()) / total_rate
            self.t += dt
            if self.p.use_nutrient:
                for _ in range(max(1, int(math.ceil(dt / self.p.nutrient_dt)))):
                    self.lattice.nutrient_step(dt)
            ev = random.choices(events, weights=rates, k=1)[0]
            self.handle_event(ev)
            self.event_count += 1
            self.record_state()
            if self.next_snapshot_idx < len(self.p.snapshot_times) and self.t >= self.p.snapshot_times[self.next_snapshot_idx]:
                self.capture_frame()
                self.next_snapshot_idx += 1

        self.capture_frame()
        end_wall = pytime.time()
        print(f"Simulation finished: t={self.t:.2f}, events={self.event_count}, wall_time={end_wall - start_wall:.1f}s")
        self.save_data()

    def handle_event(self, ev):
        etype, x, y = ev
        ct, div = self.lattice.cell_type, self.lattice.pc_div_left
        if etype == 'death':
            ct[x, y] = self.p.EMPTY; div[x, y] = 0
        elif etype == 'move':
            nb = self.lattice.get_empty_neighbors(x, y)
            if nb:
                nx, ny = random.choice(nb)
                ct[nx, ny] = ct[x, y]; div[nx, ny] = div[x, y]
                ct[x, y] = self.p.EMPTY; div[x, y] = 0
        elif etype == 'div':
            nb = self.lattice.get_empty_neighbors(x, y)
            if not nb: return
            nx, ny = random.choice(nb)
            if ct[x, y] == self.p.CSC:
                r = np.random.rand()
                if r < self.p.p_sym_self:
                    ct[nx, ny] = self.p.CSC; div[nx, ny] = 0
                elif r < self.p.p_sym_self + self.p.p_sym_diff:
                    ct[x, y] = ct[nx, ny] = self.p.PC
                    div[x, y] = div[nx, ny] = self.p.pc_max_divisions - 1
                else:
                    ct[nx, ny] = self.p.PC; div[nx, ny] = self.p.pc_max_divisions - 1
            elif ct[x, y] == self.p.PC and div[x, y] > 0:
                ct[nx, ny] = self.p.PC
                div[nx, ny] = div[x, y] - 1
                div[x, y] -= 1

    def save_data(self):
        os.makedirs(self.p.out_dir, exist_ok=True)
        np.savetxt(self.p.out_csv, self.data_records, delimiter=",",
                   header="time,CSC_count,PC_count,mean_nutrient", comments="")
        print(f"Saved data to {self.p.out_csv}")

    def save_gif(self, duration=3.0):
        os.makedirs(self.p.out_dir, exist_ok=True)
        images = [imageio.imread(f) for f in self.frames]
        imageio.mimsave(self.p.out_gif, images, duration=duration)
        print(f"GIF saved to {self.p.out_gif}")


# ======================================================
# === Command-line interface ===========================
# ======================================================

def parse_args():
    parser = argparse.ArgumentParser(description="Tumor growth simulation.")
    parser.add_argument('--config', type=str, help="Path to .ini configuration file.")
    parser.add_argument('--Lx', type=int, help="Lattice size X")
    parser.add_argument('--Ly', type=int, help="Lattice size Y")
    parser.add_argument('--t_max', type=float, help="Maximum simulation time")
    parser.add_argument('--seed', type=int, help="Random seed")
    parser.add_argument('--save_gif', action='store_true', help="Save GIF after simulation")
    return parser.parse_args()


def main():
    args = parse_args()
    p = Parameters()
    if args.config:
        p.load_from_ini(args.config)
    p.load_from_args(args)

    sim = TumorSimulation(p)
    sim.run()

    if args.save_gif:
        sim.save_gif()


if __name__ == "__main__":
    main()
