# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 11:24:20 2015

@author: rebecca
"""

from .problem import Problem
from .proxoperators import ProxOperator
from .algorithms import Algorithm
#from multiprocessing import Pool

import matplotlib, matplotlib.pyplot as pyplot

import math, scipy.fftpack, scipy.io, scipy.misc, numpy

__all__ = ["Inpainting"]



class P_new_mask(ProxOperator):
    """
    Fourier magnitude update for the Thibault et al. algorithm
    
    .. seealso:: :class:`ProxOperator`
    """
    
    def __init__(self, config):
        self.Nx = config['Nx']; self.Ny = config['Ny'];
        self.dim = math.sqrt(self.Nx*self.Ny);
        #self.ptychography_prox = config['ptychography_prox'];
        #self.fmask = config['fmask'];
        self.mask = config.mask;
        self.recon = config.recon;
        
        self.alpha = config.alpha;
    
    
    def work(self,mask):
        #Nx = self.Nx; Ny = self.Ny;
        dim = self.dim;
        c = self.mask;
        u = self.recon;
        
        # vector u-u_0+Lu
        deriv = u-self.u_0.recon + self.Lap*u;
        
        alpha = self.mask_stepsize;
        alpha = alpha * 2 * normest(spdiags(deriv*deriv,0,dim,dim))
        
        
        return Pty(phi,obj.copy(),probe.copy());



class P_new_recon(ProxOperator):
    """
    Object update for the PHeBIE algorithm
    
    .. seealso:: :class:`ProxOperator`
    """
    
    def __init__(self,config):
        self.Nx = config['Nx']; self.Ny = config['Ny'];
        self.N_pie = config['N_pie'];
        self.positions = config['positions'];
        self.overrelax = config['overrelax'];
        self.switch_object_support_constraint = config['switch_object_support_constraint'];
        self.object_support = config['object_support'];
        self.trans_max_true = config['trans_max_true'];
        self.trans_min_true = config['trans_min_true'];
    
    
    def work(self,u):
        Nx = self.Nx; Ny = self.Ny;
        positions = self.positions;
        
        rangeNx = numpy.arange(Nx,dtype=numpy.int);
        rangeNy = numpy.arange(Ny,dtype=numpy.int);
        
        obj = u.obj.copy();
        probe = u.probe;
        phi = u.phi;
        
        xsum = numpy.zeros(obj.shape);
        phisum = numpy.zeros_like(obj);
        conj_probe = numpy.conj(probe);
        modsq_probe = numpy.real(probe * conj_probe);
        
        for xpos in range(self.N_pie):
            indy = rangeNy + positions[xpos,0];
            indx = rangeNx + positions[xpos,1];
            for y in range(Ny):
                xsum[indy[y],indx] += modsq_probe[y,:];
                phisum[indy[y],indx] += conj_probe[y,:] * phi[y,:,xpos];
        
        #xsum_denom = self.overrelax * numpy.maximum(xsum,1e-30);
        xsum_denom = self.overrelax * numpy.amax(numpy.maximum( \
                                                numpy.amax(xsum,axis=0),1e-30));
        obj -= (xsum * obj - phisum) / xsum_denom;
        
        if self.switch_object_support_constraint == True:
            obj *= self.object_support;
        
        max_true = self.trans_max_true;
        min_true = self.trans_min_true;
        abs_obj = numpy.abs(obj);
        high = (abs_obj > max_true).astype(abs_obj.dtype);
        low = (abs_obj < min_true).astype(abs_obj.dtype);
        obj *= (1.0-low)*(1.0-high) + (low*min_true+high*max_true)/(abs_obj+1e-30);
        
        return Pty(phi.copy(),obj,probe.copy());



class PtychographyStats():
    """
    The Ptychography algorithms have special statistics. The necessary routines
    are collected in this class.
    """
    
    def __init__(self,config):
        self.Nx = config['Nx']; self.Ny = config['Ny'];
        self.dim = config['dim'];
        self.N_pie = config['N_pie'];
        self.norm_data = config['norm_data'];
        self.positions = config['positions'];
        self.fnorm = math.sqrt(self.Nx*self.Ny);
        self.ptychography_prox = config['ptychography_prox'];
        if not self.fmask is None:
            self.fmask = self.fmask.reshape((self.fmask.shape[0],
                                             self.fmask.shape[1],
                                             self.fmask.shape[2]*
                                             self.fmask.shape[3]));
    
    def objective(self,x1):
        """
        Get the value of the objective function (used in PALM)
        
        :param:
        x1 : Pty
            Input
        
        :return:
        number
            Objective function value of x1
        """
        obj_value = 0.0;
        phi = x1.phi;
        obj = x1.obj;
        probe = x1.probe;
        positions = self.positions;
        
        rangeNy = numpy.arange(self.Ny,dtype=numpy.int);
        rangeNx = numpy.arange(self.Nx,dtype=numpy.int);
        
        for pos in range(self.N_pie):
            indy = rangeNy + positions[pos,0];
            indx = rangeNx + positions[pos,1];
            obj_value += scipy.linalg.norm(probe * obj[indy,:][:,indx] - \
                                                        phi[:,:,pos],'fro')**2;
        return obj_value;
    
    
    def change(self, x1, x2):
        """
        Computes the change for x1 and x2
        
        :param:
        x1 : Pty
            Input 1
        x2 : Pty
            Input 2
        
        :return:
        number
            The change
        """
        norm_data = self.norm_data;
        
        changephi = 0.0;
        
        if self.Nx == 1 or self.Ny == 1:
            changephi = (scipy.linalg.norm(x1.phi-x2.phi,'fro')/norm_data)**2;
        else:
            for j in range(self.dim):
                changephi += (scipy.linalg.norm(x1.phi[:,:,j]-x2.phi[:,:,j], \
                                                            'fro')/norm_data)**2;
        changephi += (scipy.linalg.norm(x1.obj-x2.obj,'fro')/norm_data)**2;
        changephi += (scipy.linalg.norm(x1.probe-x2.probe,'fro')/norm_data)**2;
        return math.sqrt(changephi);
    
    
    def customerror(self, x1, x2, x3=None):
        """
        Computes the custom error for x1, x2 and x3
        
        Parameters
        ----------
        x1 : Pty
            Input 1
        x2 : Pty
            Input 2
        x3 : Pty, optional
            Input 3
        
        Returns
        -------
        number
            The error
        """
        object_support = self.object_support;
        sample_plane = self.sample_plane;
        positions = self.positions;
        fnorm = self.fnorm;
        fmask = self.fmask;
        amp_exp_norm = self.amp_exp_norm;
        norm_data = self.norm_data;
        
        rangeNy = numpy.arange(self.Ny,dtype=numpy.int);
        rangeNx = numpy.arange(self.Nx,dtype=numpy.int);
        
        supported_plane = sample_plane * object_support;
        
        rms_error_probe, diffphase, row_shift, col_shift = self._dftregistration( \
                scipy.fftpack.fft2(self.probe), scipy.fftpack.fft2(x2.probe));
        
        shifted_object = numpy.roll(numpy.roll(x2.obj,row_shift,axis=0), \
                                    col_shift,axis=1);
        supported_object = shifted_object * object_support;
        gammaobj = numpy.sum(supported_plane * numpy.conj(supported_object)) / \
                    numpy.sum(shifted_object*numpy.conj(shifted_object)).real;
        rms_error_obj = numpy.sum(numpy.abs(supported_plane - gammaobj * \
                        supported_object)**2) / numpy.sum(supported_plane * \
                        numpy.conj(supported_plane)).real;
                        
        r_factor = 0.0;
        for pos in range(self.N_pie):
            indy = rangeNy + positions[pos,0];
            indx = rangeNx + positions[pos,1];
            if fmask is None:
                r_factor += numpy.sum(numpy.abs(amp_exp_norm[:,:,pos] - \
                            numpy.abs(scipy.fftpack.fft2(x2.obj[indy,:][:,indx]* \
                            x2.probe)/fnorm)));
            else:
                r_factor += numpy.sum(numpy.abs(amp_exp_norm[:,:,pos] - \
                            numpy.abs(scipy.fftpack.fft2(x2.obj[indy,:][:,indx] * \
                            x2.probe)/fnorm))*fmask[:,:,pos]);
        r_factor /= numpy.sum(amp_exp_norm);
        
        change_phi = 0.0;
        if self.Nx == 1 or self.Ny == 1:
            change_phi = (scipy.linalg.norm(x1.phi-x2.phi,'fro')/norm_data)**2;
        else:
            for j in range(self.dim):
                change_phi += (scipy.linalg.norm(x1.phi[:,:,j]-x2.phi[:,:,j],'fro')/norm_data)**2;
        
        if self.ptychography_prox == 'Thibault':
            x2 = x2.copy();
            x2.phi = x3.phi.copy();
        
        objective = self.objective(x2);
        
        return (rms_error_obj, rms_error_probe, change_phi, r_factor, objective);



class PALM(Algorithm):
    """
    Proximal alternating (linearized) minimization algorithm
    
    See Also
    --------
        Algorithm : Generic interface for ProxToolbox algorithms
    """
    
    def __init__(self, config):
        """
        Parameters
        ----------
        Dictionary containing the problem configuration. It must contain the
        following mappings:
            projector_sequence:sequence of Proxoperator
                 Sequence of prox operators to be used. (The classes, no instances)
            Nx:int
                Row-dim of the product space elements
            Ny:int
                Column-dim of the product space elements
            Nz:int
                Depth-dim of the product space elements
            dim:int
                Size of the product space
            norm_data:number
                ?
            ignore_error:boolean
                Whether to ignore errors
        """
        self.projs = [p(config) for p in config['projector_sequence']];
        
        self.ignore_error = config['ignore_error'];
        self.custom_error = PtychographyStats(config);
    
    
    def run(self, u, tol, maxiter):
        """
        Runs the algorithm on some input data
        
        Parameters
        ----------
        u : array_like
            Input data
        tol : number
            Tolerance
        maxiter : int
            Maximum number of iterations
        
        Returns
        -------
        u : array_like
            Result
        u_final : array_like
            Result
        iters : int
            Number of iterations the algorithm performed
        change : array_like
            The percentage change in the norm
        gap : array_like
            Squared gap distance normalized by the magnitude constraint
        """
        projs = self.projs;
        ignore_error = self.ignore_error;
        
        iters = 0;
        change = numpy.zeros(maxiter+1);
        change[0] = 999;
        
        num_custom_errors = len(self.custom_error.customerror(u,u))
        custom_errors = numpy.zeros((num_custom_errors,maxiter+1));
        custom_errors[:,0] = 999;
        
        tmp_u = u.copy();
        
        while iters < maxiter and (ignore_error or change[iters] >= tol):
            iters += 1;
            
            for p in projs:
                tmp_u = p.work(tmp_u);
            
            if ignore_error == False:
                change[iters] = self.custom_error.change(u,tmp_u);
                custom_errors[:,iters] = self.custom_error.customerror(u,tmp_u);
                
                u = tmp_u.copy();
        
        change = change[1:iters+1];
        custom_errors = custom_errors[:,1:iters+1];
        
        return tmp_u, tmp_u.copy(), iters, change, custom_errors;



class Inpainting(Problem):
    """
    Parameters
    ----------
    """
    
    default_config = {
        'sim_data_type':'scarf.png',
        'data_path':'inputdata/',
        'Nx':100,
        'Ny':100,
        'Nz':1,
        'warmup':False,
        'warmup_alg':PALM,
        'warmup_maxiter':1,
        'algorithm':PALM,
        'maxiter':2, # This is usually 300. I reduced it to be able to test it.
        'tol':0.625,
        'ignore_error':False,
        #'ptychography_prox':'Rodenburg',
        #'blocking_switch':True,
        #'blocking_scheme':'divide',
        #'between_blocks_scheme':'sequential',
        #'within_blocks_scheme':'parallel',
        #'block_rows':2,
        #'block_cols':2,
        #'block_maxiter':1,
        #'rodenburg_inner_it':1,
        'mask_reg':0.01,
        'mask_stepsize':2.5,
        'recon_stepsize':1.2,
        'mask_inertia':0,
        'recon_inertia':0
    };
    
    
    def __init__(self,config=default_config):
        """
        Parameters
        ----------
        config : dict, optional
            Dictionary containing the problem configuration. If unspecified,
            Inpainting.default_config is used.
        """
        self.config = config.copy();
        self.config['projector_sequence'] = [P_new_mask,
                                                 P_new_recon];
        self.config['warmup_projsequence'] = [P_new_mask,
                                                  P_new_recon];
        
        
        #self.config['N_pie'] = self.config['nx']*self.config['ny'];
        self.config['dim'] = self.config['nx']*self.config['ny'];
        self.config['norm_data'] = 1;
        
    
    
    def _get_ptychography_data(self):
        """
        Generates an experiment on-the-fly. If you want to load existing
        data or make your own experiments, you can overload this.
        
        Returns
        -------
        d1x : number
            Pixel dimension
        d1y : number
            Pixel dimension
        positions : array_like
            ...
        probe : array_like
            ...
        sample_plane : array_like
            ...
        i_exp : array_like
            ...
        """
        Nx = self.config['Nx'];
        Ny = self.config['Ny'];
        nx = self.config['nx'];
        ny = self.config['ny'];
        
        rangeNx = numpy.arange(Nx,dtype=numpy.int);
        rangeNy = numpy.arange(Ny,dtype=numpy.int);
        
        # 
        # Physical parameters of the experiment
        #
        
        E = 6.2;                # X-ray energy [keV]
        l = 1e-10*12.3984/E;    # Wavelength [m]
        
        # Pixel dimensions [m]
        d2x = 172e-6;
        d2y = 172e-6;
        
        # Origin of the coordinate system
        origin_x0 = Nx/2 + 1;
        origin_y0 = Ny/2 + 1;
        
        # Position(s) of sample along optical axis
        z01 = 0.5e-4;       # distance focus <-> sample [m]
        z02 = 7.2;          # distance focus <-> detector [m]
        z12 = z02-z01;      # distance sample <-> detector [m]
        
        d1x = l*z12/(Nx*d2x);
        d1y = l*z12/(Ny*d2y);
        
        #
        # Scan parameters
        #
        
        scan_stepsize = self.config['scan_stepsize'];
        pie_step_x = scan_stepsize;
        pie_step_y = scan_stepsize;
        
        positions = None;
        
        if(self.config['scan_type'] == 'round_roi'):
            raise NotImplementedError('NYI')
        else: # i.e., scna_type == 'raster'
            positions = numpy.empty((nx*ny,2));
            n_fast = numpy.arange(nx);
            ind = 0;
            for n_slow in range(ny):
                ind = n_slow*ny;
                positions[ind:ind+nx,0] = (pie_step_y/d1y) * n_slow;
                positions[ind:ind+nx,1] = (pie_step_x/d1x) * n_fast;
        
        positions[:,0] -= numpy.amin(positions[:,0]);
        positions[:,1] -= numpy.amin(positions[:,1]);
        positions[numpy.modf(positions)[0] == 0.5] += 0.01; # numpy.round is weird
        positions = numpy.round(positions).astype(numpy.int);
        
        sample_plane = None;
        if self.config['sim_data_type'] == 'gaenseliesel':
            g = pyplot.imread('inputdata/gaenseliesel.png').astype(numpy.float);
            g /= numpy.amax(g);
            g = 0.8*g + 0.1;
            sample_plane = numpy.abs(g) * (numpy.cos(g) + 1j*numpy.sin(g));
        else:
            raise NotImplementedError('Only \'gaenseliesel\' is supported')
        
        #
        # Illumination functions
        #
        
        # FWHM of gaussian illumination [m]
        fwhm = 0.5e-6;
        w0 = fwhm/(2*math.sqrt(math.log(2)));
        
        # Get coordinate system in sample plane
        x = (numpy.arange(1,Nx+1) - ((Nx//2)+1)) * d1x;
        y = (numpy.arange(1,Ny+1) - ((Ny//2)+1)) * d1y;
        X,Y = numpy.meshgrid(x,y);
        
        # Get illumination in sample plane
        u = self._gaussu(numpy.sqrt(X**2 + Y**2),z01,l,w0);
        
        # Scan along xyz-directions
        i_exp = numpy.zeros((Ny,Nx,ny,nx));
        probe = u;
        
        # Generate the data
        ind = 0;
        for akk in range(ny):
            for bkk in range(nx):
                indy = rangeNy + positions[ind,0];
                indx = rangeNx + positions[ind,1];
                # Extract roi from sample plane
                sample = sample_plane[indy,:][:,indx];
                # Calculate exit surface wave
                esw = probe * sample;
                # Calculate coherent diffraction pattern
                i_exp[:,:,akk,bkk] = numpy.abs(numpy.rot90(numpy.roll( \
                    numpy.roll(scipy.fftpack.fft2(esw),origin_x0-1,axis=0), \
                    origin_y0-1,axis=1),-2))**2;
                
                ind += 1;
        
        self.config['trans_min_true'] = 0;
        self.config['trans_max_true'] = 1;
        
        return d1x, d1y, l, z01, positions, probe, sample_plane, i_exp;
    
    
    def _create_blocks(self):
        """
        Create blocks for the current problem
        """
        sample_plane = self.config['sample_plane'];
        probe = self.config['probe'];
        positions = self.config['positions'];
        
        rows = self.config['block_rows'];
        cols = self.config['block_cols'];
        
        Nx = self.config['Nx'];
        Ny = self.config['Ny'];
        N_pie = self.config['N_pie'];
        
        rangeNy = numpy.arange(Ny);
        rangeNx = numpy.arange(Nx);
        
        # Use the initial probe to estimate the object support
        probe_support = (numpy.real(numpy.conj(probe) * probe) > 0.01).astype(numpy.float);
        shifted_probe_supports = numpy.zeros(sample_plane.shape+tuple([N_pie]));
        objectsupport = numpy.zeros_like(sample_plane);
        for pos in range(N_pie):
            indy = rangeNy + positions[pos,0];
            indx = rangeNx + positions[pos,1];
            
            for y in range(Ny):
                shifted_probe_supports[indy[y],indx,pos] = probe_support[y,:];
                objectsupport[indy[y],indx] += probe_support[y,:];
            objectsupport = (objectsupport > 0.01).astype(objectsupport.dtype);
        #probes = shifted_probe_supports;
        
        # Apply the selected blocking strategy
        placed = numpy.zeros(N_pie,dtype=numpy.bool);
        in_block = numpy.zeros(N_pie,dtype=numpy.int);
        blocking_scheme = self.config['blocking_scheme'];
        if blocking_scheme == 'one':
            in_block.fill(0);
        elif blocking_scheme == 'single_view':
            placed.fill(True);
            in_block = numpy.arange(N_pie);
        elif blocking_scheme == 'divide':
            the_max = numpy.amax(positions,axis=0)+1;
            row_len = the_max[0] / rows;
            col_len = the_max[1] / cols;
            for p in range(N_pie):
                top_index = positions[p,:];
                for r in range(rows):
                    for c in range(cols):
                        b = c + r*cols;
                        if placed[p] == False and \
                           top_index[0] >= round(r*row_len) and \
                           top_index[0] <= round((r+1)*row_len) and \
                           top_index[1] >= round(c*col_len) and \
                           top_index[1] <= round((c+1)*col_len):
                            in_block[p] = b;
                            placed[p] = True;
        elif blocking_scheme == 'distribute':
            in_block = numpy.arange(N_pie) % (rows*cols);
        else:
            raise NotImplementedError('NYI')
        self.config['in_block'] = in_block;
    
    
    def _presolve(self):
        Nx = self.config['Nx'];
        Ny = self.config['Ny'];
        nx = self.config['nx'];
        ny = self.config['ny'];
        
        i_exp = None;
        N_pie = self.config['N_pie'];
        
        rangeNx = numpy.arange(Nx,dtype=numpy.int);
        rangeNy = numpy.arange(Ny,dtype=numpy.int);
        
        d1x, d1y, l, z01, positions, probe, sample_plane, i_exp =  \
            self._get_ptychography_data();
        
        #
        # Process the data (based on Pär's data processor)
        #
        
        # Average total intensity in experimental diffraction pattern
        i_sum = numpy.sum(numpy.sum(i_exp,axis=1),axis=0);
        #i_avg = numpy.mean(numpy.mean(i_sum,axis=0));
        #print 'Average intensity of single frame is %3.2e photons' % i_avg;
        # Maximum average intensity in a single experimental frame
        i_max_avg = numpy.amax(i_sum) / (Nx*Ny);
        
        # Generate cell array of normalized intensities
        # Normalization such, that for I = I_max (as matrix) the average pixel
        # intensity is 1 and total is Nx*Ny, for all other I the value is lower
        amp_exp_norm = numpy.zeros((i_exp.shape[0],i_exp.shape[1],N_pie))
        
        # Specimen plane (E1) / probe plane (as long as direct propagation is used)
        X,Y = numpy.meshgrid((rangeNx-(Nx//2)-1)*d1x,(rangeNy-(Ny//2)-1)*d1y);
        
        # Relative radius of ideal circular probe mask with respect to radius of
        # pinhole used for initial simulation of the probe function
        r = 2.2e-6;
        beta_sam = 0.9*Nx*d1x/2/r * self.config['probe_mask_gamma'];
        
        probe_mask = None;
        if self.config['switch_probemask'] == True:
            probe_mask = (numpy.sign(-(Y**2 + X**2 - (beta_sam*r)**2)) + 1.0)/2;
        else:
            probe_mask = numpy.ones(probe.shape);
        
        # Generate noise
        if not self.config['noise'] is None:
            raise NotImplementedError('NYI')
        
        cx = i_exp.shape[1]//2;
        cy = i_exp.shape[0]//2;
        bs_mask = self.config['bs_mask'];
        bs_factor = self.config['bs_factor'];
        fmask = self.config['fmask'];
        # Removal of the residual noise in the data is done here
        ind = 0;
        for n_slow in range(ny):
            for n_fast in range(nx):
                # Check if semi-transparent central stop was used. Scale
                # intensity accordingly.
                if not bs_mask is None:
                    dat = i_exp[:,:,n_slow,n_fast];
                    dat[bs_mask] *= bs_factor;
                    i_exp[:,:,n_slow,n_fast] = dat;
                
                # Save time and get experimental data into MATLAB's non-centered
                # system.
                # I wonder if this is necessary with Scipy
                dat = numpy.rot90(numpy.roll(numpy.roll( \
                    i_exp[:,:,n_slow,n_fast],-cy,axis=0),-cx,axis=1),2);
                
                # Normalization such that average intensity of data is on the
                # order of 1
                amp_exp_norm[:,:,ind] = numpy.sqrt(dat/i_max_avg);
                
                # Check for pixelmask
                if not fmask is None:
                    fmask[:,:,n_slow,n_fast] = numpy.rot90(numpy.roll(numpy.roll( \
                        fmask[:,:,n_slow,n_fast],-cy,axis=0),-cx,axis=1),2);
                
                ind += 1;
        
        guess_value = 1000;
        
        # Initial guess for the probe
        probe_guess_type = self.config['probe_guess_type'];
        probe_guess = None;
        if probe_guess_type == 'exact':
            probe_guess = probe.copy();
        elif probe_guess_type == 'exact_amp':
            probe_guess = numpy.abs(probe);
        elif probe_guess_type == 'circle':
            probe_guess = numpy.zeros_like(probe);
            radius = 13.0/256.0 * probe.shape[0];
            x_c = probe.shape[0]/2;
            y_c = probe.shape[1]/2;
            for i in range(probe.shape[0]):
                for j in range(probe.shape[1]):
                    if math.sqrt((i-x_c)**2 + (j-y_c)**2) < radius:
                        probe_guess[i,j] = guess_value;
        elif probe_guess_type == 'robin_initial':
            X,Y = numpy.meshgrid((rangeNx-(Nx//2)-1)*d1x,(rangeNy-(Ny//2)-1)*d1y);
            def heaviside(M):
                return (1+numpy.sign(M))/2;
            def ellipsis(X,Y,x0,y0,a,b):
                return heaviside(-(((X-x0)/a)**2+((Y-y0)/b)**2-1));
            probe_guess = numpy.abs(scipy.misc.imrotate(ellipsis(X,Y,0,0,0.2e-6,0.2e-6),0));
            probe_guess = self._prop_nf_pa(probe_guess,l,z01,d1x,d1y);
            probe_guess /= scipy.linalg.norm(probe_guess,'fro') /           \
                                math.sqrt(numpy.size(probe_guess));
        elif probe_guess_type == 'ones':
            probe_guess = numpy.ones_like(probe);
        else:
            raise NotImplementedError('NYI')
        
        # Initial guess for the object
        object_guess_type = self.config['object_guess_type'];
        object_guess = None;
        if object_guess_type == 'exact':
            object_guess = sample_plane.copy();
        elif object_guess_type == 'ones':
            object_guess = numpy.ones_like(sample_plane) / guess_value;
        elif object_guess_type == 'constant':
            object_guess = numpy.ones_like(sample_plane) / guess_value;
        elif object_guess_type == 'exact_perturbed':
            object_guess = sample_plane + 2e-1 * \
                            (numpy.random.random_sample(sample_plane.shape)-0.5);
        elif object_guess_type == 'random':
            object_guess = (1.0/guess_value) * \
                            numpy.random.random_sample(sample_plane.shape) * \
                            numpy.exp(1j*(2 * numpy.pi * \
                                numpy.random.random_sample(sample_plane.shape) - \
                                numpy.pi)/2);
        else:
            raise NotImplementedError('NYI')
        
        
        # Compute object support constraint (using initial probe guess)
        low = 0.1;
        mask = None;
        object_support = numpy.zeros(sample_plane.shape);
        if probe_guess_type == 'robinGaussian':
            mask = probe_mask;
        else:
            mask = ((numpy.abs(probe_guess)**2)>low).astype(object_support.dtype);
        
        for pos in range(N_pie):
            indy = rangeNy + positions[pos,0];
            indx = rangeNx + positions[pos,1];
            for y in range(Ny):
                object_support[indy[y],indx] += mask[y,:];
        
        object_support = (object_support > 0).astype(object_support.dtype);
        
        object_guess *= object_support;
        sample_plane *= object_support;
        
        if self.config['switch_probemask'] == True:
            probe *= probe_mask;
        
        # Generate initial guess for exit wave functions
        phi = numpy.ones((Nx,Ny,N_pie),dtype=probe_guess.dtype);
        for pos in range(N_pie):
            indy = rangeNy + positions[pos,0];
            indx = rangeNx + positions[pos,1];
            phi[:,:,pos] = probe_guess * object_guess[indy,:][:,indx];
        
        self.u = Pty(phi,object_guess,probe_guess);
        
        self.config['probe'] = probe;
        self.config['sample_plane'] = sample_plane;
        self.config['amp_exp_norm'] = amp_exp_norm;
        self.config['probe_mask'] = probe_mask;
        self.config['positions'] = positions;
        self.config['object_support'] = object_support;
        self.config['cfact'] = 0.0001*nx*ny;
        
        
        #
        # Generate blocks for the Ptychography problem using a specified strategy
        #
        if self.config['blocking_switch'] == True:
            self._create_blocks();
        
        
        # Put together our algorithms
        if self.config['blocking_switch'] == True:
            self.algorithm = Blocking_meta(self.config);
        else:
            self.algorithm = self.config['algorithm'](self.config);
        
        if self.config['warmup'] == True:
            self.warmup_config = self.config.copy();
            self.warmup_config['algorithm'] = self.config['warmup_alg'];
            self.warmup_config['projector_sequence'] = self.config['warmup_projsequence'];
            if self.config['blocking_switch'] == True:
                self.warmup_alg = Blocking_meta(self.warmup_config);
            else:
                self.warmup_alg = self.config['warmup_alg'](self.warmup_config);
        


    def _solve(self):
        if self.config['warmup'] == True:
            print ('Starting warm-up algorithm')
            u, self.u_warm, iters, gap, change = \
                self.warmup_alg.run(self.u,self.config['tol'],self.config['warmup_maxiter']);
            print ('End warm-up algorithm')
        else:
            print ('Skipping warm-up algorithm')
            self.u_warm = self.u;
        
        print ('Start algorithm')
        u, self.u_final, self.iters, self.change, self.custom_errors = \
                self.algorithm.run(self.u_warm,self.config['tol'],self.config['maxiter']);
        print ('End algorithm')
    
    
    def _postsolve(self):
        obj = self.u_final.obj;
        probe = self.u_final.probe;
        
        sample_plane = self.config['sample_plane'];
        positions = self.config['positions'];
        probe_mask = self.config['probe_mask'];
        object_support = self.config['object_support'];
        N_pie = self.config['N_pie'];
        Nx = self.config['Nx'];
        Ny = self.config['Ny'];
        
        rangeNx = numpy.arange(Nx,dtype=numpy.int);
        rangeNy = numpy.arange(Ny,dtype=numpy.int);
        
        pyplot.figure('Object amplitude');
        ax = pyplot.subplot(1,2,1);
        ax.title.set_text('Best approximation');
        pyplot.imshow(numpy.abs(obj));
        ax = pyplot.subplot(1,2,2);
        ax.title.set_text('Exact object');
        pyplot.imshow(numpy.abs(sample_plane));
        
        pyplot.figure('Object phase');
        ax = pyplot.subplot(1,2,1);
        ax.title.set_text('Best approximation');
        pyplot.imshow(numpy.angle(obj));
        ax = pyplot.subplot(1,2,2);
        ax.title.set_text('Exact object');
        pyplot.imshow(numpy.angle(sample_plane));
        
        # Probe
        pyplot.figure('Probe');
        ax = pyplot.subplot(1,2,1);
        ax.title.set_text('Best probe approximation');
        pyplot.imshow(matplotlib.colors.hsv_to_rgb(self._im2hsv(probe,1.0)));
        ax = pyplot.subplot(1,2,2);
        ax.title.set_text('Exact probe');
        pyplot.imshow(matplotlib.colors.hsv_to_rgb(self._im2hsv(self.config['probe'],1.0)));
        
        pyplot.figure('Objective value');
        pyplot.semilogy(self.custom_errors[4,:]);
        pyplot.xlabel('Iteration');
        pyplot.ylabel('Objective function value');
        
        meas = numpy.zeros(sample_plane.shape);
        for pos in range(N_pie):
            indy = rangeNy + positions[pos,0];
            indx = rangeNx + positions[pos,1];
            for y in range(Ny):
                meas[indy[y],indx] += probe_mask[y,:];
        pyplot.figure('Number of measurements, pixelwise');
        pyplot.imshow(meas);
        
        pyplot.figure('Change');
        pyplot.plot(self.change);
        pyplot.xlabel('Iteration');
        pyplot.ylabel('change');
        
        pyplot.figure('RMS Error Object');
        ax = pyplot.subplot(1,2,1);
        ax.title.set_text('RMS Error Object');
        ax.xaxis.set_label('Iteration');
        ax.yaxis.set_label('RMS-Error, Object');
        pyplot.semilogy(self.custom_errors[0,:]);
        ax = pyplot.subplot(1,2,2);
        ax.title.set_text('Area over which error is measured');
        pyplot.imshow(object_support);
        
        pyplot.figure('RMS Error Probe');
        pyplot.semilogy(self.custom_errors[1,:]);
        pyplot.xlabel('Iteration');
        pyplot.ylabel('RMS-Error, Probe');
        
        pyplot.figure('Change in Phi');
        pyplot.semilogy(self.custom_errors[2,:]);
        pyplot.xlabel('Iteration');
        pyplot.ylabel('Change in Phi');
        
        pyplot.figure('R-Factor');
        pyplot.plot(self.custom_errors[3,:]);
        pyplot.xlabel('Iteration');
        pyplot.ylabel('R-Factor');
        
        pyplot.show();

