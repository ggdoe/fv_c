import numpy as np
import matplotlib.pyplot as plt

# Aliases for the variables
ddl = 4
IR = 0
IMU = 1
IU = 1
IMV = 2
IV = 2
IP = 3
IE = 3

# Main class
class Sim:
    def __init__(self, 
                problem='SOD', 
                N=(128, 128), 
                xy_max=(1.0,1.0),
                Ng=2, 
                BC='absorbing', 
                riemann='hll', 
                muscl=False,
                CFL=0.8, 
                recons='plm', 
                plot_all=False, 
                wave_estimate='davis',
                gamma=1.4, 
                limiter='minmod',
                tend=1.0,
                max_ite=1000, 
                put_obstacle=False,
                obstacle_xy=(32, 40, 32, 40)):
        '''
        Constructor of the Sim class

        Parameters :
        ------------
            problem (string): type of problem to solve
            Nx (int): number of cells
            Ng (int): number of ghosts
            BC (string): type of boundary conditions
            riemann (string): type of riemann solver
            CFL (float): Courant condition factor (0<CFL<1)
            recons (string): type of reconstruction
            plot_all (boolean): should the system plot each iteration of the code
            wave_estimate (string) : which method should be used for the wave-speed estimates
            gamma (float): adiabatic index
            limiter (string): type of slope-limiter to apply
            tend (float): maximum time of the simulation
            max_ite (int): maximum number of iterations 
        '''
        self.problem = problem
        self.CFL = CFL
        self.recons = recons
        self.limiter = limiter
        self.gamma = gamma
        self.riemann = riemann
        self.BC = BC
        self.plot_all = plot_all
        self.wave_estimate = wave_estimate
        self.tend = tend
        self.max_ite = max_ite
        self.muscl = muscl

        # Creating the grid
        self.Nx, self.Ny = N
        self.Ng = Ng
        self.Ntot_x = self.Nx + 2*Ng
        self.Ntot_y = self.Ny + 2*Ng
        xmax, ymax = xy_max
        self.dx = xmax / self.Nx
        self.dy = ymax / self.Ny
        self.lo = Ng
        self.hix = Ng+self.Nx-1
        self.hiy = Ng+self.Ny-1

        self.x = np.linspace((0.5 - Ng)*self.dx, xmax +
                             (Ng-0.5)*self.dx, self.Ntot_x)
        self.y = np.linspace((0.5 - Ng)*self.dy, ymax +
                             (Ng-0.5)*self.dy, self.Ntot_y)
        self.U = np.zeros((ddl, self.Ntot_y, self.Ntot_x))
        self.Q = np.zeros((ddl, self.Ntot_y, self.Ntot_x))
        self.slopes_x = np.zeros((ddl, self.Ntot_y, self.Ntot_x))
        self.slopes_y = np.zeros((ddl, self.Ntot_y, self.Ntot_x))

        self.xdom = self.x[self.lo:self.hix+1]
        self.ydom = self.y[self.lo:self.hiy+1]
        self.Udom = self.U[:, self.lo:self.hiy+1, self.lo:self.hix+1]
        self.Qdom = self.Q[:, self.lo:self.hiy+1, self.lo:self.hix+1]

        self.init_problem()

        if put_obstacle:
            xmin, xmax, ymin, ymax = obstacle_xy
            self.Qdom[:, ymin:ymax+1, xmin:xmax+1] = 0.0
        self.Q0 = np.copy(self.Qdom)

        self.ite        = 0
        self.t          = 0.0
        self.dt         = 1.0e-2
        self.max_dt_var = 1.1

        self.put_obstacle = put_obstacle
        self.obstacle_xy = obstacle_xy
        self.tmp = False

    def init_problem(self):
        '''
        Initializing domain with a few common problems
        '''
        if self.problem == 'SOD':
            side = len(self.x) // 2
            self.Q[IR, :side, :] = 1.0
            self.Q[IR, side:, :] = 0.125

            self.Q[IP, :side, :] = 1.0
            self.Q[IP, side:, :] = 0.1

            # shape = np.shape(self.Q[IU, :, :])
            # self.Q[IU, :, :] = 0.5 * (2*np.random.rand(*shape) - 1)
            # self.Q[IV, :, :] = 0.5 * (2*np.random.rand(*shape) - 1)

        elif self.problem == 'SOD-2d':
            side = len(self.x) // 2
            self.Q[IR, :side, :side] = 1.0
            self.Q[IR, side:, :side] = 0.125*2
            self.Q[IR, side:, side:] = 1.0*1.5
            self.Q[IR, :side, side:] = 0.125

            self.Q[IP, :side, :side] = 1.0
            self.Q[IP, side:, :side] = 0.1/1.5
            self.Q[IP, side:, side:] = 1.0/2
            self.Q[IP, :side, side:] = 0.1

        elif self.problem == 'laxliu3':
            # voir -> Lax, Liu 1998
            ix = len(self.x) // 2
            iy = len(self.y) // 2
            Q_ = (self.Q[:, :iy, :ix], self.Q[:, :iy, ix:], self.Q[:, iy:, :ix], self.Q[:, iy:, ix:])
            
            r_ = (1.5, 0.5323, 0.138, 0.5323)
            u_ = (0, 1.206, 1.206, 0)
            v_ = (0, 0, 1.206, 1.206)
            p_ = (1.5, 0.3, 0.029, 0.3)

            for Q, r, u, v, p in zip(Q_, r_, u_, v_, p_):
                Q[IR,:] = r
                Q[IU,:] = u
                Q[IV,:] = -v
                Q[IP,:] = p

        elif self.problem == 'kelvin-helmholtz':
            ly = np.fabs(self.y - 0.5) > 0.25
            ry = ~ly
            self.Q[IR, ly, :] = 2.0
            self.Q[IR, ry, :] = 1.0
            self.Q[IU, ly, :] = -0.5
            self.Q[IU, ry, :] = 0.5
            self.Q[IV,  :, :] = 0.0
            self.Q[IP,  :, :] = 2.5

            # np.random.seed(1)
            # shape = np.shape(self.Q[IU, :, :])
            # self.Q[IU, :, :] += 0.01 * (np.random.rand(*shape) - 0.5)
            # self.Q[IV, :, :] += 0.01 * (np.random.rand(*shape) - 0.5)
            
            lx = self.y < 0.5
            self.Q[IU, :,  lx] += +0.01
            self.Q[IV, :, ~lx] += -0.01

        elif self.problem == 'gaussian':
            r2 = np.matrix((self.x - 0.5)**2)
            self.Q[IR] = np.exp(-100.0*(r2.T + r2)) + 1.0
            self.Q[IU]   = +1.0
            self.Q[IV]   = 0.5
            self.Q[IP]   = 1.0

        elif self.problem == 'implosion':
            x = np.matrix(self.x)
            y = np.matrix(self.y)
            id = (x + y.T > 0.5)
            self.Q[IR, id] = 1.0
            self.Q[IP, id] = 1.0
            self.Q[IR, ~id] = 0.125
            self.Q[IP, ~id] = 0.14
            self.Q[IU,  :, :] = 0.0
            self.Q[IV,  :, :] = 0.0

        # elif self.problem == 'stationary_wave':
        #     left = self.x < 0.5
        #     right = ~left

        #     self.Q[IR, left] = 1.4
        #     self.Q[IR, right] = 1.0
        #     self.Q[IU, :] = 0.0
        #     self.Q[IP, :] = 1.0


        self.primitive_to_conservative()
        self.iteration = 0

    def primitive_to_conservative(self):
        ''' Converts from Q to U '''
        self.U[IR] = self.Q[IR]
        self.U[IMU] = self.Q[IR] * self.Q[IU]
        self.U[IMV] = self.Q[IR] * self.Q[IV]

        Ek = 0.5 * self.Q[IR] * (self.Q[IU]**2.0 + self.Q[IV]**2.0)
        self.U[IE] = Ek + self.Q[IP] / (self.gamma-1.0)

    def conservative_to_primitive(self):
        ''' Converts from U to Q'''
        self.Q[IR] = self.U[IR]
        self.Q[IU] = self.U[IMU] / self.U[IR]
        self.Q[IV] = self.U[IMV] / self.U[IR]
        Ek = 0.5 * (self.U[IMU]**2.0 + self.U[IMV]**2.0) / self.U[IR]
        self.Q[IP] = (self.U[IE] - Ek) * (self.gamma-1.0)

        # Checking that everything is alright
        neg_dens = self.Q[IR] < 0.0
        neg_pres = self.Q[IP] < 0.0

        if np.sum(neg_dens) + np.sum(neg_pres) > 0:
            print("neg!!!")
            self.Q[IR, neg_dens] = 1.0e-5
            self.Q[IP, neg_pres] = 1.0e-5
            self.primitive_to_conservative()

    def get_primitive_to_conservative(self, Q):
            U = np.zeros_like(Q)
            U[IR] = Q[IR]
            U[IMU] = Q[IR] * Q[IU]
            U[IMV] = Q[IR] * Q[IV]

            Ek = 0.5 * Q[IR] * (Q[IU]**2.0 + Q[IV]**2.0)
            U[IE] = Ek + Q[IP] / (self.gamma-1.0)
            return U

    def get_conservative_to_primitive(self, U):
            Q = np.zeros_like(U)
            Q[IR] = U[IR]
            Q[IU] = U[IMU] / U[IR]
            Q[IV] = U[IMV] / U[IR]

            Ek = 0.5 * Q[IR] * (Q[IU]**2.0 + Q[IV]**2.0)
            Q[IP] = (U[IE] - Ek) * (self.gamma - 1)
            return Q

    def fill_boundaries(self):
        '''
        Fills the boundaries according to what is asked from the simulation
        '''
        Ng = self.Ng; lo = self.lo; hix = self.hix; hiy = self.hiy

        if self.BC == 'absorbing':
            # top/bot/left/right
            self.U[:, :lo,    lo:hix+1] = np.tile(self.U[:,lo,lo:hix+1],Ng).reshape(ddl, Ng, hix+1-lo)
            self.U[:, hiy+1:, lo:hix+1] = np.tile(self.U[:,hiy,lo:hix+1],Ng).reshape(ddl, Ng, hix+1-lo)
            self.U[:, lo:hiy+1,    :lo] = np.repeat(self.U[:,lo:hiy+1,lo],Ng).reshape(ddl, hiy+1-lo, Ng)
            self.U[:, lo:hiy+1, hix+1:] = np.repeat(self.U[:,lo:hiy+1,hix],Ng).reshape(ddl, hiy+1-lo, Ng)
            # TODO : transverse

        elif self.BC == 'reflect':
            '''by repeat last row/col'''
            # self.U[:, :lo,    lo:hix+1] = np.tile(self.U[:,lo,lo:hix+1],Ng).reshape(ddl, Ng, hix+1-lo)
            # self.U[:, hiy+1:, lo:hix+1] = np.tile(self.U[:,hiy,lo:hix+1],Ng).reshape(ddl, Ng, hix+1-lo)
            # self.U[:, lo:hiy+1,    :lo] = np.repeat(self.U[:,lo:hiy+1,lo],Ng).reshape(ddl, hiy+1-lo, Ng)
            # self.U[:, lo:hiy+1, hix+1:] = np.repeat(self.U[:,lo:hiy+1,hix],Ng).reshape(ddl, hiy+1-lo, Ng)
            '''by reflecting lasts rows/cols'''
            self.U[:, :lo,    lo:hix+1] = self.U[:, lo+Ng-1:lo-1:-1, lo:hix+1]
            self.U[:, hiy+1:, lo:hix+1] = self.U[:, hiy:hiy-Ng:-1,   lo:hix+1]
            self.U[:, lo:hiy+1,    :lo] = self.U[:, lo:hiy+1, lo+Ng-1:lo-1:-1]
            self.U[:, lo:hiy+1, hix+1:] = self.U[:, lo:hiy+1,   hix:hix-Ng:-1]

            # reflect orthogonal speed
            self.U[IMU, lo:hiy+1,    :lo] = -self.U[IMU, lo:hiy+1,    :lo]   
            self.U[IMU, lo:hiy+1, hix+1:] = -self.U[IMU, lo:hiy+1, hix+1:] 
            self.U[IMV, :lo,    lo:hix+1] = -self.U[IMV, :lo,    lo:hix+1]   
            self.U[IMV, hiy+1:, lo:hix+1] = -self.U[IMV, hiy+1:, lo:hix+1] 

        elif self.BC == 'periodic':
            self.U[:,     :lo , lo:hix+1] = self.U[:, hiy-Ng+1:hiy+1,      lo:hix+1]   # left
            self.U[:,   hiy+1:, lo:hix+1] = self.U[:,       lo:lo+Ng,      lo:hix+1]   # right
            self.U[:, lo:hiy+1,      :lo] = self.U[:,      lo:hiy+1, hix-Ng+1:hix+1]   # top
            self.U[:, lo:hiy+1,   hix+1:] = self.U[:,        lo:hiy+1,     lo:lo+Ng]   # bottom

            self.U[:,     :lo,      :lo] = self.U[:, hiy-Ng+1:hiy+1, hix-Ng+1:hix+1]   # top-left
            self.U[:,   hiy+1:,     :lo] = self.U[:,     lo:lo+Ng, hix-Ng+1:hix+1]   # top-right
            self.U[:,   hiy+1:,  hix+1:] = self.U[:,     lo:lo+Ng,     lo:lo+Ng]   # bottom-right
            self.U[:,     :lo,   hix+1:] = self.U[:, hiy-Ng+1:hiy+1,     lo:lo+Ng]   # bottom-left

    def obstacle(self, xmin, xmax, ymin, ymax):
        ''' put rectangular obstacle in the conservative vector U'''
        Ng = self.Ng
        rect = self.Udom[:, ymin-Ng:ymax+1+Ng, xmin-Ng:xmax+1+Ng]
        lo = Ng; hiy = ymax - ymin + 1 + Ng; hix = xmax - xmin + 1 + Ng

        rect[:, lo:lo+Ng,       lo:hix+1] = rect[:, :lo,    lo:hix+1]
        rect[:, hiy+1-Ng:hiy+1, lo:hix+1] = rect[:, hiy+1:, lo:hix+1]
        rect[:, lo:hiy+1, lo:lo+Ng      ] = rect[:, lo:hiy+1,    :lo]
        rect[:, lo:hiy+1, hix+1-Ng:hix+1] = rect[:, lo:hiy+1, hix+1:]

        # reflect orthogonal speed
        # rect[IMU, lo:lo+Ng,       lo:hix+1] = -rect[IMU, lo:lo+Ng,       lo:hix+1]
        # rect[IMU, hiy+1-Ng:hiy+1, lo:hix+1] = -rect[IMU, hiy+1-Ng:hiy+1, lo:hix+1]
        # rect[IMV, lo:hiy+1, lo:lo+Ng      ] = -rect[IMV, lo:hiy+1,       lo:lo+Ng]
        # rect[IMV, lo:hiy+1, hix+1-Ng:hix+1] = -rect[IMV, lo:hiy+1, hix+1-Ng:hix+1]

        # remove orthogonal speed
        rect[IMU, lo:lo+Ng,       lo:hix+1] = 0.0
        rect[IMU, hiy+1-Ng:hiy+1, lo:hix+1] = 0.0
        rect[IMV, lo:hiy+1, lo:lo+Ng      ] = 0.0
        rect[IMV, lo:hiy+1, hix+1-Ng:hix+1] = 0.0

    def compute_dt(self):
        ''' 
        Computes the time-step according to the CFL condition
        '''
        cs = np.sqrt(self.Q[IP]*self.gamma / self.Q[IR])

        # new_dt = self.CFL * np.min(self.dx / (cs + np.fabs(self.Q[IU]) + np.fabs(self.Q[IV])))
        min_x = np.min(self.dx / (cs + np.fabs(self.Q[IU])))
        min_y = np.min(self.dy / (cs + np.fabs(self.Q[IV])))
        new_dt = self.CFL * min(min_x, min_y)

        # norm_v = np.sqrt(self.Q[IU]**2 + self.Q[IV]**2)
        # new_dt = self.CFL * np.min(self.dx / (cs + norm_v))


        # Limiting the maximum value of dt to avoid too big variations
        # if new_dt/self.dt > self.max_dt_var:
        #     new_dt = self.dt * self.max_dt_var
        
        # Limiting dt if we go over the max sim time
        # new_dt = min(self.tend-self.t+1.0e-6, new_dt)
        self.dt = new_dt

    def compute_slopes(self):
        '''
        Computes the slope according to a given slope limiter
        '''
        def minmod(a, b):
            select_a = np.fabs(a) < np.fabs(b)
            select_b = ~select_a
            res = np.zeros_like(a)
            res[select_a] = a[select_a]
            res[select_b] = b[select_b]
            res[a*b < 0.0] = 0.0
            return res
        
        self.slopes_x[:,:,:] = 0.0
        self.slopes_y[:,:,:] = 0.0
        for i in range(ddl):
            if self.limiter == 'nolimiter':
                self.slopes_x[i, :, 1:-1] = 0.5 * (self.Q[i, :, 2:] - self.Q[i, :, :-2])
                self.slopes_y[i, 1:-1, :] = 0.5 * (self.Q[i, 2:, :] - self.Q[i, :-2, :])

            elif self.limiter == 'minmod':
                dl = self.Q[i, :, 1:-1] - self.Q[i, :,  :-2]
                dr = self.Q[i, :,   2:] - self.Q[i, :, 1:-1]
                self.slopes_x[i, :, 1:-1] = minmod(dl, dr)

                dl = self.Q[i, 1:-1, :] - self.Q[i, :-2, :]
                dr = self.Q[i,   2:, :] - self.Q[i, 1:-1, :]
                self.slopes_y[i, 1:-1, :] = minmod(dl, dr)

            elif self.limiter == 'monotonized_central':
                dl = self.Q[i, :, 1:-1] - self.Q[i, :,  :-2]
                dr = self.Q[i, :,   2:] - self.Q[i, :, 1:-1]
                self.slopes_x[i, :, 1:-1] = minmod(0.5*(dl+dr), 2.0*minmod(dl, dr))

                dl = self.Q[i, 1:-1, :] - self.Q[i,  :-2, :]
                dr = self.Q[i,   2:, :] - self.Q[i, 1:-1, :]
                self.slopes_y[i, 1:-1, :] = minmod(0.5*(dl+dr), 2.0*minmod(dl, dr))

    def get_flux_from_p(self, Q):
        ''' Computes Euler fluxes of a given primitive array'''
        rho = Q[IR]
        p   = Q[IP]
        u   = Q[IU]
        v   = Q[IV]

        Ek = 0.5 * rho * (u**2.0 + v**2.0)
        E = Ek + p / (self.gamma-1.0)

        fluxes = np.zeros_like(Q)
        fluxes[IR] = rho * u
        fluxes[IU] = rho * u * u + p
        fluxes[IV] = rho * u * v
        fluxes[IE] = (E + p) * u

        return fluxes

    def compute_wave_speed_estimates(self, QL, QR):
        aL = np.sqrt(QL[IP] * self.gamma / QL[IR])
        aR = np.sqrt(QR[IP] * self.gamma / QR[IR])

        if self.wave_estimate == 'davis':
            sl_minL = QL[IU] - aL
            sl_maxL = QL[IU] + aL
            sl_minR = QR[IU] - aR
            sl_maxR = QR[IU] + aR

            SL = np.minimum(sl_minL, sl_minR)
            SR = np.maximum(sl_maxL, sl_maxR)


        elif self.wave_estimate == 'toro2019':
            z = (self.gamma-1.0) / (2.0*self.gamma)
            ps = ((aL + aR - self.gamma*z*(QR[IU]-QL[IU])) / (aL/QL[IP]**z + aR/QR[IP]**z))**(1.0/z)

            z = (self.gamma+1.0) / (2.0*self.gamma)
            qL = np.where(ps <= QL[IP], 1.0, np.sqrt(1.0 + z * (ps/QL[IP]-1.0)))
            qR = np.where(ps <= QR[IP], 1.0, np.sqrt(1.0 + z * (ps/QR[IP]-1.0)))
            
            SL = QL[IU] - aL*qL
            SR = QR[IU] + aR*qR

        return SL, SR

    def solve_fluxes(self, QL, QR):
        ''' 
        Applies the riemann solver with left values QL and right values QR
        The returned array fluxes is of shape (ddl, Ntot-2). fluxes[:, i] are the fluxes for interface i+1/2
        '''

        fluxes = np.zeros_like(QL)
        if self.riemann == 'hll':
            SL, SR = self.compute_wave_speed_estimates(QL, QR)

            # Inefficient because we compute the fluxes regardless of the signal speed 
            # estimates
            # TODO : improve that.
            FL = self.get_flux_from_p(QL)
            FR = self.get_flux_from_p(QR)

            apply_left = 0.0 <= SL
            apply_hll = (SL <= 0.0) & (0.0 <= SR)
            apply_right = SR <= 0.0

            # Rebuilding conservative values for HLL flux
            # Note : This could be optimized by only writing using the apply_hll mask
            UL = self.get_primitive_to_conservative(QL)
            UR = self.get_primitive_to_conservative(QR)
            
            FHLL = (SR*FL - SL*FR + SL*SR*(UR-UL)) / (SR-SL)

            fluxes[:, apply_left]  = FL[:, apply_left]
            fluxes[:, apply_right] = FR[:, apply_right]
            fluxes[:, apply_hll]   = FHLL[:, apply_hll]

        elif self.riemann == 'hllc':
            SL, SR = self.compute_wave_speed_estimates(QL, QR)
            UL = self.get_primitive_to_conservative(QL)
            UR = self.get_primitive_to_conservative(QR)

            Sstar = (QR[IP] - QL[IP] + \
                     UL[IMU] * (SL - QL[IU]) - UR[IMU] * (SR - QR[IU])) / \
                     (QL[IR] * (SL - QL[IU]) - QR[IR] * (SR - QR[IU]))
            
            def Ustar(U, S, Q):
                Us = np.zeros_like(U)
                Us[IR,:] = 1.0
                Us[IMU,:] = Sstar
                Us[IMV,:] = Q[IV]
                Us[IE,:] = U[IE] / U[IR] + (Sstar - Q[IU]) * (Sstar + Q[IP] / (Q[IR] * (S - Q[IU])))
                Z = Q[IR] * (S - Q[IU]) / (S - Sstar)
                return Us * Z
            
            FL = self.get_flux_from_p(QL)
            FR = self.get_flux_from_p(QR)
            FLstar = FL + SL * (Ustar(UL, SL, QL) - UL)
            FRstar = FR + SR * (Ustar(UR, SR, QR) - UR)

            apply_left       = 0.0 <= SL
            apply_left_star  = (SL <= 0.0) & (0.0 <= Sstar)
            apply_right_star = (Sstar <= 0.0) & (0.0 <= SR)
            apply_right      = SR <= 0.0

            fluxes[:, apply_left]       = FL[:, apply_left]
            fluxes[:, apply_left_star]  = FLstar[:, apply_left_star]
            fluxes[:, apply_right_star] = FRstar[:, apply_right_star]
            fluxes[:, apply_right]      = FR[:, apply_right]

        return fluxes
    
    def muscl_update(self, q):
        dx = self.slopes_x
        dy = self.slopes_y
        dt = self.dt
        out = np.zeros_like(q)

        dtdx = 0.5 * dt / self.dx
        dtdy = 0.5 * dt / self.dy

        out[IR] = q[IR] + (- q[IU] * dx[IR] - q[IR] * dx[IU]) * dtdx \
                        + (- q[IV] * dy[IR] - q[IR] * dy[IV]) * dtdy
        out[IU] = q[IU] + (- q[IU] * dx[IU] - dx[IP] / q[IR]) * dtdx \
                        + (- q[IV] * dy[IU]) * dtdy
        out[IV] = q[IV] + (- q[IU] * dx[IV]) * dtdx \
                        + (- q[IV] * dy[IV] - dy[IP] / q[IR]) * dtdy
        out[IP] = q[IP] + (- self.gamma * q[IP] * dx[IU] - q[IU] * dx[IP]) * dtdx \
                        + (- self.gamma * q[IP] * dy[IV] - q[IV] * dy[IP]) * dtdy
        
        return out[:,:,:]

    def update_cells(self):
        '''
        Does the update of each cell in the domain
        '''

        # 0. half time update
        if self.muscl:
            Q = self.muscl_update(self.Q)
        else:
            Q = self.Q

        # 1. 
        '''Reconstructing'''
        if self.recons == 'plm':
            xL = Q[:, :, :-1] + 0.5 * self.slopes_x[:, :, :-1]  # xL[i] -> Left interface i+1/2
            xR = Q[:, :,  1:] - 0.5 * self.slopes_x[:, :,  1:]  # xR[i] -> Right interface i+1/2
            yL = Q[:, :-1, :] + 0.5 * self.slopes_y[:, :-1, :]  # yL[i] -> Left interface i+1/2
            yR = Q[:,  1:, :] - 0.5 * self.slopes_y[:,  1:, :]  # yR[i] -> Right interface i+1/2
        elif self.recons == 'pcm':
            xL = np.copy(Q[:, :, :-1])
            xR = np.copy(Q[:, :,  1:])
            yL = np.copy(Q[:, :-1, :])
            yR = np.copy(Q[:,  1:, :])
        # 2.

        # swap u,v direction for the riemann solver
        yL[IU], yL[IV] = yL[IV], np.copy(yL[IU])
        yR[IU], yR[IV] = yR[IV], np.copy(yR[IU])

        '''Getting fluxes x'''
        F = self.solve_fluxes(xL, xR)
        # print(F[IU])
        FL = F[:, :, :-1]
        FR = F[:, :,  1:]

        '''Getting fluxes y'''
        G = self.solve_fluxes(yL, yR)
        G[IU], G[IV] = G[IV], np.copy(G[IU])
        GL = G[:, :-1, :]
        GR = G[:,  1:, :]
        # print(G[IU])

        # 3. Updating cells
        self.U[:, :, 1:-1] += self.dt / self.dx * (FL-FR)
        self.U[:, 1:-1:, ] += self.dt / self.dy * (GL-GR)

    def step(self):
        '''
        One step of the Godunov Solver
        '''

        if self.put_obstacle:
            self.obstacle(*self.obstacle_xy)


        # 3. Computing dt using the CFL condition
        self.compute_dt()
        print(self.dt)

        # 4. Calculating slopes
        self.compute_slopes()

        # 5. Solving riemann problems and updating cells
        self.update_cells()
        # 1. Filling boundaries
        self.fill_boundaries()
        # 2. Convert to primitives
        self.conservative_to_primitive()

    def plot(self):
        '''
        Plots the data to a .png file
        '''
        self.conservative_to_primitive()
        plt.close('all')
        fig, ax = plt.subplots(1, 3, figsize=(15, 7))
        ax[0].plot(self.xdom, self.Q0[IR], '--k')
        ax[1].plot(self.xdom, self.Q0[IU], '--k')
        ax[2].plot(self.xdom, self.Q0[IP], '--k')
        ax[0].plot(self.xdom, self.Qdom[IR], '-+r')
        ax[1].plot(self.xdom, self.Qdom[IU], '-+r')
        ax[2].plot(self.xdom, self.Qdom[IP], '-+r')

        plt.savefig('plt.{:06}.png'.format(self.ite))

    def run(self, field_to_plot=IR, vlim=(0, 0), cmap='viridis', freq_save = 1/60.):
        '''
        Running the simulation up to a maximum of iterations or a certain time
        '''
        self.t = 0.0
        last_save = 0.0
        if vlim[0] == vlim[1]:
            vmin, vmax = np.min(self.Qdom[field_to_plot]), np.max(self.Qdom[field_to_plot])
            d = vmax - vmin
            vmin = vmin - 0.20 * d
            vmax = vmax + 0.20 * d
        else:
            vmin, vmax = vlim[0], vlim[1]
        # plt.imsave(f'out/tmp/{0:06}.png', self.Q0[field_to_plot], vmin=vmin, vmax=vmax, cmap=cmap)


        np.set_printoptions(precision=3, linewidth=20000, threshold=20000, formatter={'float': '{:0.3f}'.format})
        # print(self.Q0[IR])
        for i in range(1000):
            self.ite += 1
            self.step();
            plt.imsave(f'out/{self.ite:06}.png', self.Q[IR], vmin=0.8, vmax=2.2, cmap=cmap)
        # self.step();
        # self.step();
        # self.step();self.conservative_to_primitive()
        # print(self.dt)
        print(self.Q[IR])
        # print(self.slopes_x[IR])
        # plt.imshow(self.Qdom[IR])
        # plt.imshow(self.Q0[IR])
        plt.show()
        return;

        while self.t < self.tend:
            self.step()
            self.ite += 1
            self.t += self.dt
            if self.plot_all:
                self.plot()
            
            # if self.t > last_save + freq_save :
            #     last_save = self.t
            #     self.conservative_to_primitive()
                
            #     if self.put_obstacle:
            #         xmin, xmax, ymin, ymax = self.obstacle_xy
            #         self.Qdom[:, ymin:ymax+1, xmin:xmax+1] = 0.0
                # plt.imsave(f'out/tmp/{self.ite:06}.png', self.Qdom[field_to_plot], vmin=vmin, vmax=vmax, cmap=cmap)
                # plt.ylim((0.1,1.025))
                # plt.plot(self.ydom, self.Qdom[field_to_plot, self.Nx // 2, :])
                # plt.savefig(f'out/tmp/{self.ite:06}.png')
                # plt.close()

            # print(f'Iteration {self.ite}\tt={self.t:.5}\tdt={self.dt:.5}\t{"*" if last_save == self.t else ""}')
            
            if self.ite == self.max_ite:
                break
        print(self.Q[IR])

if __name__ == '__main__':
    import os
    
    """ param """
    n = 64
    N = (n,n)
    xy_max=(1.0,1.0)
    # problem='gaussian'
    problem='kelvin-helmholtz'
    # problem='implosion'
    # problem='SOD'
    # problem='SOD-2d' 
    # BC='absorbing'
    BC='periodic' 
    # BC='reflect' 
    # recons='pcm' 
    recons='plm'
    limiter='minmod'
    # limiter='nolimiter'
    # limiter='monotonized_central'
    CFL=0.3
    max_ite=0
    tend=0.2
    wave_estimate='toro2019'
    riemann = 'hllc'
    # wave_estimate='davis'
    # riemann = 'hll'

    put_obstacle=False
    # put_obstacle=False
    obstacle_xy=(32, 40, 32, 40)
    # obstacle_xy=(64, 80, 64, 80)

    """ plot """
    field_to_plot=IR
    # vlim = (0.8, 2.2)
    vlim = (0.2, 1.5)
    # vlim = (0,0)
    cmap='nipy_spectral_r'
    freq_save = 1/60.
    muscl=False
    
    sim = Sim(N=N, xy_max=xy_max, problem=problem, BC=BC, recons=recons, limiter=limiter, CFL=CFL, max_ite=max_ite, tend=tend, wave_estimate=wave_estimate, riemann=riemann, put_obstacle=put_obstacle, obstacle_xy=obstacle_xy, muscl=muscl)
    # os.system(f"rm out/tmp/*.png")
    sim.run(field_to_plot=field_to_plot, vlim=vlim, cmap=cmap, freq_save=freq_save)
    # os.system(f"ffmpeg -y -pattern_type glob -i 'out/tmp/*.png' -vf 'scale=640:-1' -c:v libx264 out/out.mp4")
    
    # plt.plot(sim.xdom, sim.Qdom[IR, sim.Nx // 2, :])
    # plt.show()