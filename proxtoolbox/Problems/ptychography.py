# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 15:20:03 2014

@author: stefan

'ptychography' is an instance of the Problems-module.
It implements all tools required for solving the Ptychography-problem
described in the following paper

    References
    ----------
    .. [1] R. Hesse, D.R. Luke, S. Sabach, and M.K. Tam, "Proximal
       Heterogeneous Block Implicit-Explicit Method and Application to Blind
       Ptychographic Diffraction Imaging", SIAM J. Imaging Sciences,
       8(1):426–457, 2015.
"""

from .problems import Problem
# from proxtoolbox.Algorithms import PALM
from proxtoolbox.Algorithms.algorithms import Algorithm
from proxtoolbox.ProxOperators import magproj
from proxtoolbox.ProxOperators.proxoperators import ProxOperator
from multiprocessing import Pool

from matplotlib import colors, pyplot

from  math import atan, atan2, ceil, exp, log, sqrt
import scipy.fftpack, scipy.io, scipy.misc
import numpy

#__all__ = ["Ptychography","Ptychography_NTT_01_26210"]



class Pty:
    """
    In the Ptychography problem phi, the object and the probe do not have the
    same dimensions. This class combines them into a single data structure.
    Some standard operators are overloaded for convenience.
    
    """
    
    def __init__(self,phi,obj,probe):
        """
        Initialization
        
        Parameters
        ----------
        phi : array_like
              The p-th diffraction pattern
        obj : array_like
              The object
        probe : array_like
                The probe
        """
        self.phi = phi;
        self.obj = obj;
        self.probe = probe;
        
        self.dtype = self.phi.dtype;
    
    def __add__(self,p):
        """
        Elementwise addition
        """
        return Pty(self.phi+p.phi, self.obj+p.pbj, self.probe+p.probe);
    
    def __sub__(self,p):
        """
        Elementwise subtraction
        """
        return Pty(self.phi-p.phi, self.obj-p.pbj, self.probe-p.probe);
    
    def __rmul__(self,t):
        """
        Multiplication
        """
        return Pty(self.phi*t, self.obj*t, self.probe*t);
    
    def __div__(self,t):
        """
        Division
        """
        return Pty(self.phi/t, self.obj, self.probe);
    
    def __getitem__(self,key):
        """
        Getter for an item of phi
        """
        return self.phi[key];
    
    def copy(self):
        """
        Creates a copy of a Pty-instance
        """
        return Pty(self.phi.copy(),self.obj.copy(),self.probe.copy())
    
    



class P_rodenburg(ProxOperator):
    """
    Rodenburg's Ptychography algorithm
    
    Notes
    -----
    Done by
        * Pär Mattsson, Universität Göttingen, 2013
        * Matthew Tam, CARMA Centre, University of Newcastle, 2014
    """
    
    def __init__(self,config):
        """
        Initialization
        """
        self.Nx = config['Nx']; self.Ny = config['Ny'];
        self.N_pie = config['N_pie'];
        self.positions = config['positions'];
        self.switch_probemask = config['switch_probemask'];
        self.probemask = config['probe_mask'];
        self.amp_exp_norm = config['amp_exp_norm'];
        self.inner_it = config['rodenburg_inner_it'];
        self.switch_object_support_constraint = config['switch_object_support_constraint'];
        self.object_support = config['object_support'];
        self.trans_max_true = config['trans_max_true'];
        self.trans_min_true = config['trans_min_true'];
        self.fmask = config['fmask'];
        if not self.fmask is None:
            self.fmask = self.fmask.reshape((self.fmask.shape[0],
                                             self.fmask.shape[1],
                                             self.fmask.shape[2]*
                                             self.fmask.shape[3]));
        self.fnorm = sqrt(self.Nx*self.Ny);
    
        
    def work(self, u):
        """
        Parameters
        ----------
        u : Pty - Ptychography has special needs, so u must be a Pty.
        
        Returns
        -------
        Pty - The output iterate
        
        See Also
        --------
        :class:`.Pty`
        """
        fnorm = self.fnorm;
        positions = self.positions;
        Nx = self.Nx; Ny = self.Ny;
        switch_probemask = self.switch_probemask;
        probemask = self.probemask;
        amp_exp_norm = self.amp_exp_norm;
        inner_it = self.inner_it;
        switch_object_support_constraint = self.switch_object_support_constraint;
        object_support = self.object_support;
        trans_max_true = self.trans_max_true;
        trans_min_true = self.trans_min_true;
        fmask = self.fmask;
        
        this_object = u.obj.copy();
        this_probe = u.probe.copy();
        this_phi = numpy.zeros_like(u.phi);
        
        rangeNx = numpy.arange(Nx,dtype=numpy.int);
        rangeNy = numpy.arange(Ny,dtype=numpy.int);
        
        for pos in numpy.random.permutation(self.N_pie):
            indy = rangeNy + positions[pos,0];
            indx = rangeNx + positions[pos,1];
            
            if switch_probemask == True:
                probe_norm = scipy.linalg.norm(this_probe,'fro') / \
                                sqrt(numpy.size(this_probe));
                this_probe *= probemask;
                this_probe *= probe_norm/(scipy.linalg.norm(this_probe,'fro') / \
                                sqrt(numpy.size(this_probe)));
            
            phi = this_probe * this_object[indy,:][:,indx];
            old_phi_hat = scipy.fftpack.fft2(phi) / fnorm;
            phi_hat = magproj(old_phi_hat,amp_exp_norm[:,:,pos]);
            if not fmask is None:
                phi_hat = phi_hat * fmask[:,:,pos] + old_phi_hat * \
                            (fmask[:,:,pos] == 0).astype(fmask.dtype);
            phi_prime = scipy.fftpack.ifft2(phi_hat) * fnorm;
            
            for k in range(inner_it):
                this_object_old = this_object.copy();
                this_probe_old = this_probe.copy();
                cthis_probe_old = numpy.conj(this_probe_old);
                cthis_object_old = numpy.conj(this_object_old);
                i_probe = numpy.amax(cthis_probe_old*this_probe_old) + 1e-8;
                i_object = numpy.amax(this_object_old[indy,:][:,indx] * \
                            cthis_object_old[indy,:][:,indx]) + 1e-8;
                d_phi = phi_prime - phi;
                
                n_probe = (cthis_object_old[indy,:][:,indx]/i_object) * d_phi;
                n_object = (cthis_probe_old/i_probe) * d_phi;
                
                this_probe = this_probe_old + n_probe;
                for y in range(Ny):
                    this_object[indy[y],indx] = \
                                this_object_old[indy[y],indx] + n_object[y,:];
                
                phi = this_probe * this_object[indy,:][:,indx];
                phi_hat = scipy.fftpack.fft2(phi) / fnorm;
                phi_hat = magproj(phi_hat,amp_exp_norm[:,:,pos]);
                phi_prime = scipy.fftpack.ifft2(phi_hat) * fnorm;
                
            if switch_object_support_constraint == True:
                this_object *= object_support;
            
            abs_object = numpy.abs(this_object);
            high = (abs_object > trans_max_true).astype(this_phi.dtype);
            low = (abs_object < trans_min_true).astype(this_phi.dtype);
            this_object = ((1-low)*(1-high)*this_object) + \
                ((low*trans_min_true)+(high*trans_max_true)) * \
                this_object / (abs_object+1e-30);
            
            this_phi[:,:,pos] = this_probe * this_object[indy,:][:,indx];
        
        return Pty(this_phi,this_object,this_probe);


class P_rodenburg_probe_fixed(P_rodenburg):
    """
    Rodenburg's Ptychography algorithm
    
    """

    def work(self,u):
        """
        Parameters
        ----------
        u : Pty - Input
        """
        fnorm = self.fnorm;
        positions = self.positions;
        Nx = self.Nx; Ny = self.Ny;
        switch_probemask = self.switch_probemask;
        probemask = self.probemask;
        amp_exp_norm = self.amp_exp_norm;
        inner_it = self.inner_it;
        switch_object_support_constraint = self.switch_object_support_constraint;
        object_support = self.object_support;
        trans_max_true = self.trans_max_true;
        trans_min_true = self.trans_min_true;
        fmask = self.fmask;
        
        this_object = u.obj.copy();
        this_probe = u.probe.copy();
        this_phi = numpy.zeros_like(u.phi);
        
        rangeNx = numpy.arange(Nx,dtype=numpy.int);
        rangeNy = numpy.arange(Ny,dtype=numpy.int);
        
        for pos in numpy.random.permutation(range(self.N_pie)):
            indy = rangeNy + positions[pos,0];
            indx = rangeNx + positions[pos,1];
            
            if switch_probemask == True:
                probe_norm = scipy.linalg.norm(this_probe,'fro') / \
                                sqrt(numpy.size(this_probe));
                this_probe *= probemask;
                this_probe *= probe_norm/(scipy.linalg.norm(this_probe,'fro') / \
                                sqrt(numpy.size(this_probe)));
            
            phi = this_probe * this_object[indy,:][:,indx];
            old_phi_hat = scipy.fftpack.fft2(phi) / fnorm;
            phi_hat = magproj(old_phi_hat,amp_exp_norm[:,:,pos]);
            if not fmask is None:
                phi_hat = phi_hat * fmask[:,:,pos] + old_phi_hat * \
                            (fmask[:,:,pos] == 0).astype(fmask.dtype);
            phi_prime = scipy.fftpack.ifft2(phi_hat) * fnorm;
            
            for k in range(inner_it):
                this_object_old = this_object.copy();
                this_probe_old = this_probe.copy();
                cthis_probe_old = numpy.conj(this_probe_old);
                i_probe = numpy.amax(numpy.real(cthis_probe_old*this_probe_old))+1e-8;
                #i_probe = numpy.amax(numpy.abs(this_probe_old)**2) + 1e-8;
                d_phi = phi_prime - phi;
                
                n_object = (cthis_probe_old/i_probe) * d_phi;
                
                for y in range(Ny):
                    this_object[indy[y],indx] = \
                                this_object_old[indy[y],indx] + n_object[y,:];
                
                phi = this_probe * this_object[indy,:][:,indx];
                phi_hat = scipy.fftpack.fft2(phi) / fnorm;
                phi_hat = magproj(phi_hat,amp_exp_norm[:,:,pos]);
                phi_prime = scipy.fftpack.ifft2(phi_hat) * fnorm;
                
            if switch_object_support_constraint == True:
                this_object *= object_support;
            
            abs_object = numpy.abs(this_object);
            high = (abs_object > trans_max_true).astype(this_phi.dtype);
            low = (abs_object < trans_min_true).astype(this_phi.dtype);
            this_object = ((1.0-low)*(1.0-high)*this_object) + \
                ((low*trans_min_true)+(high*trans_max_true)) * \
                this_object / (abs_object+1e-30);
            
            this_phi[:,:,pos] = this_probe * this_object[indy,:][:,indx];
        
        return Pty(this_phi,this_object,this_probe);


class P_thibault_f(ProxOperator):
    """
    Fourier magnitude update for the Thibault et al. algorithm
    
    """
    
    def __init__(self, config):
        """
        Initialization
        """
        self.Nx = config['Nx']; self.Ny = config['Ny'];
        self.N_pie = config['N_pie'];
        self.positions = config['positions'];
        self.amp_exp_norm = config['amp_exp_norm'];
        self.fnorm = sqrt(self.Nx*self.Ny);
        self.ptychography_prox = config['ptychography_prox'];
        self.fmask = config['fmask'];
        if not self.fmask is None:
            self.fmask = self.fmask.reshape((self.fmask.shape[0],
                                             self.fmask.shape[1],
                                             self.fmask.shape[2]*
                                             self.fmask.shape[3]));
    
    
    def work(self,u):
        """
        Parameters
        ----------
        u : Pty - Input
        """
        Nx = self.Nx; Ny = self.Ny;
        fnorm = self.fnorm;
        positions = self.positions;
        amp_exp_norm = self.amp_exp_norm;
        fmask = self.fmask;
        
        phi = u.phi.copy();
        obj = u.obj;
        probe = u.probe;
        
        if self.ptychography_prox == 'Thibault_AP':
            rangeNx = numpy.arange(Nx,dtype=numpy.int);
            rangeNy = numpy.arange(Ny,dtype=numpy.int);
            for pos in range(self.N_pie):
                indy = rangeNy + positions[pos,0];
                indx = rangeNx + positions[pos,1];
                phi[:,:,pos] = probe * obj[indy,:][:,indx];
        
        for pos in range(self.N_pie):
            old_phi_hat = scipy.fftpack.fft2(phi[:,:,pos]) / fnorm;
            phi_hat = magproj(old_phi_hat,amp_exp_norm[:,:,pos]);
            if not fmask is None:
                phi_hat = phi_hat * fmask[:,:,pos] + old_phi_hat * \
                            (fmask[:,:,pos] == 0).astype(fmask.dtype);
            phi[:,:,pos] = scipy.fftpack.ifft2(phi_hat) * fnorm;
        
        return Pty(phi,obj.copy(),probe.copy());


class P_thibault_o(ProxOperator):
    """
    The object update for the Thibault et al. algorithm, with probe fixed
    
    """
    
    def __init__(self, config):
        """
        Initialization
        """
        self.Nx = config['Nx']; self.Ny = config['Ny'];
        self.N_pie = config['N_pie'];
        self.positions = config['positions'];
        self.switch_object_support_constraint = config['switch_object_support_constraint'];
        self.object_support = config['object_support'];
        self.trans_max_true = config['trans_max_true'];
        self.trans_min_true = config['trans_min_true'];
        self.switch_probemask = config['switch_probemask'];
        self.probe_mask = config['probe_mask'];
        self.ptychography_prox = config['ptychography_prox'];
        
    
    def work(self, u):
        """
        Parameters
        ----------
        u : Pty - Input
        """
        Nx = self.Nx; Ny = self.Ny;
        N_pie = self.N_pie;
        positions = self.positions;
        switch_object_support_constraint = self.switch_object_support_constraint;
        object_support = self.object_support;
        trans_max_true = self.trans_max_true;
        trans_min_true = self.trans_min_true;
        switch_probemask = self.switch_probemask;
        probe_mask = self.probe_mask;
        
        rangeNx = numpy.arange(Nx,dtype=numpy.int);
        rangeNy = numpy.arange(Ny,dtype=numpy.int);
        
        phi = u.phi;
        
        this_probe = u.probe.copy();
        this_object = None;
        
        obj_enum = numpy.empty_like(u.obj);
        obj_denom = numpy.empty_like(u.obj);
        
        for i in (0,1,2):
            obj_enum.fill(1e-30);
            obj_denom.fill(1e-30);
            c_probe = numpy.conj(this_probe)
            i_probe = numpy.real(c_probe*this_probe);
            for pos in range(N_pie):
                indy = rangeNy + positions[pos,0];
                indx = rangeNx + positions[pos,1];
                
                for y in range(Ny):
                    obj_enum[indy[y],indx] += c_probe[y,:] * phi[y,:,pos];
                    obj_denom[indy[y],indx] += i_probe[y,:];
            
            this_object = obj_enum / obj_denom;
            
            if switch_object_support_constraint == True:
                this_object *= object_support;
            
            abs_object = numpy.abs(this_object);
            high = (abs_object > trans_max_true);
            low = (abs_object < trans_min_true);
            this_object *= (1.0-low)*(1.0-high) + \
                            (low*trans_min_true+high*trans_max_true) /        \
                            (abs_object+1e-30);
            
            if switch_probemask == True:
                this_probe *= probe_mask;
        
        phi = phi.copy();
        
        if self.ptychography_prox == 'Thibault':
            for pos in range(N_pie):
                indy = rangeNy + positions[pos,0];
                indx = rangeNx + positions[pos,1];
                phi[:,:,pos] = this_probe * this_object[indy,:][:,indx];
        
        return Pty(phi,this_object,this_probe);


class P_thibault_op(ProxOperator):
    """
    The simultaneous probe and object update for the Thibault et al. algorithm
    
    """
    
    def __init__(self,config):
        """
        Initialization
        """
        self.Nx = config['Nx']; self.Ny = config['Ny'];
        self.N_pie = config['N_pie'];
        self.positions = config['positions'];
        self.switch_object_support_constraint = config['switch_object_support_constraint'];
        self.object_support = config['object_support'];
        self.trans_max_true = config['trans_max_true'];
        self.trans_min_true = config['trans_min_true'];
        self.switch_probemask = config['switch_probemask'];
        self.probe_mask = config['probe_mask'];
        self.cfact = config['cfact'];
        self.ptychography_prox = config['ptychography_prox'];
        
    
    def work(self, u):
        """
        Parameters
        ----------
        u : Pty - Input
        """
        Nx = self.Nx; Ny = self.Ny;
        N_pie = self.N_pie;
        positions = self.positions;
        switch_object_support_constraint = self.switch_object_support_constraint;
        object_support = self.object_support;
        trans_max_true = self.trans_max_true;
        trans_min_true = self.trans_min_true;
        switch_probemask = self.switch_probemask;
        probe_mask = self.probe_mask;
        cfact = self.cfact;
        
        rangeNx = numpy.arange(Nx,dtype=numpy.int);
        rangeNy = numpy.arange(Ny,dtype=numpy.int);

        phi = u.phi;
        
        this_probe = u.probe.copy();
        this_object = None;
        
        obj_enum = numpy.empty_like(u.obj);
        obj_denom = numpy.empty_like(u.obj);
        
        for i in (0,1,2,3,4):
            obj_enum.fill(1e-8);
            obj_denom.fill(1e-8);
            c_probe = numpy.conj(this_probe);
            i_probe = numpy.real(c_probe*this_probe);
            
            for pos in range(N_pie):
                indy = rangeNy + positions[pos,0];
                indx = rangeNx + positions[pos,1];
                
                for y in range(Ny):
                    obj_enum[indy[y],indx] += c_probe[y,:] * phi[y,:,pos];
                    obj_denom[indy[y],indx] += i_probe[y,:];
                
            this_object = obj_enum / obj_denom;
            
            if switch_object_support_constraint == True:
                this_object *= object_support;
            
            abs_object = numpy.abs(this_object);
            high = (abs_object > trans_max_true);
            low = (abs_object < trans_min_true);
            this_object *= (1.0-low)*(1.0-high) + \
                        (low*trans_min_true+high*trans_max_true) /        \
                        (abs_object+1e-30);
            
            probe_enum = cfact * this_probe;
            probe_denom = cfact;
            conj_obj = numpy.conj(this_object);
            abs_obj = numpy.real(this_object*conj_obj);
            for pos in range(N_pie):
                indy = rangeNy + positions[pos,0];
                indx = rangeNx + positions[pos,1];
                
                probe_enum += conj_obj[indy,:][:,indx] * phi[:,:,pos];
                probe_denom += abs_obj[indy,:][:,indx];
            
            this_probe = probe_enum / probe_denom;
            
            if switch_probemask == True:
                probe_norm = scipy.linalg.norm(this_probe,'fro') /        \
                                sqrt(numpy.size(this_probe));
                this_probe *= probe_mask;
                this_probe *= probe_norm /      \
                                (scipy.linalg.norm(this_probe,'fro') /    \
                                sqrt(numpy.size(this_probe)));
        
        phi = phi.copy();
        
        if self.ptychography_prox == 'Thibault':
            for pos in range(N_pie):
                indy = rangeNy + positions[pos,0];
                indx = rangeNx + positions[pos,1];
                phi[:,:,pos] = this_probe * this_object[indy,:][:,indx];
            
        return Pty(phi,this_object,this_probe);


class P_PHeBIE_probe(ProxOperator):
    """
    Probe update for the PHeBIE algorithm
    
    """
    
    def __init__(self,config):
        """
        Initialization
        """
        self.Nx = config['Nx']; self.Ny = config['Ny'];
        self.N_pie = config['N_pie'];
        self.positions = config['positions'];
        self.overrelax = config['overrelax'];
        self.switch_probemask = config['switch_probemask'];
        self.probe_mask = config['probe_mask'];
        
    def work(self,u):
        """
        Parameters
        ----------
        u : Pty - Input
        """
        positions = self.positions;
        
        rangeNx = numpy.arange(self.Nx,dtype=numpy.int);
        rangeNy = numpy.arange(self.Ny,dtype=numpy.int);
        
        phi = u.phi;
        probe = u.probe.copy();
        
        ysum = numpy.zeros(probe.shape);
        phisum = numpy.zeros_like(probe);
        
        conj_obj = numpy.conj(u.obj);
        modsq_obj = numpy.real(u.obj * conj_obj);
        
        for ypos in range(self.N_pie):
            indy = rangeNy + positions[ypos,0];
            indx = rangeNx + positions[ypos,1];
            ysum += modsq_obj[indy,:][:,indx];
            phisum += conj_obj[indy,:][:,indx] * phi[:,:,ypos];
        
        ysum_denom = self.overrelax * numpy.amax(numpy.maximum( \
                                                numpy.amax(ysum,axis=0),1e-30));
        probe -= (ysum * probe - phisum) / ysum_denom;
        
        if self.switch_probemask == True:
            probe *= self.probe_mask;
        
        abs_probe = numpy.abs(probe);
        high = (abs_probe > 10e30).astype(abs_probe.dtype);
        obj = None;
        if numpy.amax(high) == 1:
            obj = (1.0-high) * u.obj + (high*10e30) * probe / (abs_probe+1e-30);
        else:
            obj = u.obj.copy();
        
        return Pty(phi.copy(), obj, probe);


class P_PHeBIE_probe_ptwise(P_PHeBIE_probe):
    """
    Probe update for the PHeBIE algorithm
    
    """
    
    def work(self,u):
        """
        Parameters
        ----------
        u : Pty - Input
        """
        positions = self.positions;
        
        rangeNx = numpy.arange(self.Nx,dtype=numpy.int);
        rangeNy = numpy.arange(self.Ny,dtype=numpy.int);
        
        phi = u.phi;
        probe = u.probe.copy();
        
        ysum = numpy.zeros(probe.shape);
        phisum = numpy.zeros_like(probe);
        
        conj_obj = numpy.conj(u.obj);
        modsq_obj = numpy.real(u.obj * conj_obj);
        
        for ypos in range(self.N_pie):
            indy = rangeNy + positions[ypos,0];
            indx = rangeNx + positions[ypos,1];
            ysum += modsq_obj[indy,:][:,indx];
            phisum += conj_obj[indy,:][:,indx] * phi[:,:,ypos];
        
        #ysum_denom = self.overrelax * numpy.amax(numpy.maximum( \
        #                                        numpy.amax(ysum,axis=0),1e-30));
        ysum_denom = self.overrelax * numpy.maximum(ysum,1e-30);
        probe -= (ysum * probe - phisum) / ysum_denom;
        
        if self.switch_probemask == True:
            probe *= self.probe_mask;
        
        abs_probe = numpy.abs(probe);
        high = (abs_probe > 10e30).astype(abs_probe.dtype);
        obj = None;
        if numpy.any(high) == True:
            obj = (1.0-high) * u.obj + (high*10e30) * probe / (abs_probe+1e-30);
        else:
            obj = u.obj.copy();
        
        return Pty(phi.copy(), obj, probe);


class P_PHeBIE_object(ProxOperator):
    """
    Object update for the PHeBIE algorithm
    
    """
    
    def __init__(self,config):
        """
        Initialization
        """
        self.Nx = config['Nx']; self.Ny = config['Ny'];
        self.N_pie = config['N_pie'];
        self.positions = config['positions'];
        self.overrelax = config['overrelax'];
        self.switch_object_support_constraint = config['switch_object_support_constraint'];
        self.object_support = config['object_support'];
        self.trans_max_true = config['trans_max_true'];
        self.trans_min_true = config['trans_min_true'];
    
    
    def work(self,u):
        """
        Parameters
        ----------
        u : Pty - Input
        """
        positions = self.positions;
        
        rangeNx = numpy.arange(self.Nx,dtype=numpy.int);
        rangeNy = numpy.arange(self.Ny,dtype=numpy.int);
        
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
            for y in range(self.Ny):
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


class P_PHeBIE_object_ptwise(P_PHeBIE_object):
    """
    Object update for the PHeBIE algorithm
    
    """
    
    def work(self,u):
        """
        Parameters
        ----------
        u : Pty - Input
        """
        positions = self.positions;
        
        rangeNx = numpy.arange(self.Nx,dtype=numpy.int);
        rangeNy = numpy.arange(self.Ny,dtype=numpy.int);
        
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
            for y in range(self.Ny):
                xsum[indy[y],indx] += modsq_probe[y,:];
                phisum[indy[y],indx] += conj_probe[y,:] * phi[y,:,xpos];
        
        xsum_denom = self.overrelax * numpy.maximum(xsum,1e-30);
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


class P_PHeBIE_phi(ProxOperator):
    """
    Phi update for the PHeBIE algorithm
    
    """
    
    def __init__(self,config):
        """
        Initialization
        """
        self.Nx = config['Nx']; self.Ny = config['Ny'];
        self.N_pie = config['N_pie'];
        self.fnorm = sqrt(self.Nx*self.Ny);
        self.positions = config['positions'];
        self.amp_exp_norm = config['amp_exp_norm'];
        self.fmask = config['fmask'];
        if not self.fmask is None:
            self.fmask = self.fmask.reshape((self.fmask.shape[0],
                                             self.fmask.shape[1],
                                             self.fmask.shape[2]*
                                             self.fmask.shape[3]));
        
    
    def work(self,u):
        """
        Parameters
        ----------
        u : Pty - Input
        """
        fnorm = self.fnorm;
        positions = self.positions;
        amp_exp_norm = self.amp_exp_norm;
        fmask = self.fmask;
        
        rangeNx = numpy.arange(self.Nx,dtype=numpy.int);
        rangeNy = numpy.arange(self.Ny,dtype=numpy.int);
        
        phi = u.phi.copy();
        obj = u.obj;
        probe = u.probe;
        
        for pos in range(self.N_pie):
            indy = rangeNy + positions[pos,0];
            indx = rangeNx + positions[pos,1];
            new_phi = probe * obj[indy,:][:,indx];
            old_phi_hat = scipy.fftpack.fft2(new_phi) / fnorm;
            phi_hat = magproj(old_phi_hat,amp_exp_norm[:,:,pos]);
            if not fmask is None:
                phi_hat = phi_hat * fmask[:,:,pos] + old_phi_hat * \
                            (fmask[:,:,pos] == 0).astype(fmask.dtype);
            phi[:,:,pos] = scipy.fftpack.ifft2(phi_hat) * fnorm;
        
        return Pty(phi,obj.copy(),probe.copy());
    

class P_PHeBIE_phi_regularized(ProxOperator):
    """
    Regularized phi update for the PHeBIE algorithm
    
    """
    
    def __init__(self,config):
        """
        Initialization
        """
        self.Nx = config['Nx']; self.Ny = config['Ny'];
        self.N_pie = config['N_pie'];
        self.fnorm = sqrt(self.Nx*self.Ny);
        self.positions = config['positions'];
        self.overrelax = config['overrelax'];
        self.amp_exp_norm = config['amp_exp_norm'];
        self.fmask = config['fmask'];
        if not self.fmask is None:
            self.fmask = self.fmask.reshape((self.fmask.shape[0],
                                             self.fmask.shape[1],
                                             self.fmask.shape[2]*
                                             self.fmask.shape[3]));
        
    
    def work(self,u):
        """
        Parameters
        ----------
        u : Pty - Input
        """
        fnorm = self.fnorm;
        positions = self.positions;
        amp_exp_norm = self.amp_exp_norm;
        fmask = self.fmask;
        
        rangeNx = numpy.arange(self.Nx,dtype=numpy.int);
        rangeNy = numpy.arange(self.Ny,dtype=numpy.int);
        
        phi = u.phi.copy();
        obj = u.obj;
        probe = u.probe;
        
        gamma = self.overrelax - 1;
        
        for pos in range(self.N_pie):
            indy = rangeNy + positions[pos,0];
            indx = rangeNx + positions[pos,1];
            new_phi = (2.0/(2.0+gamma)) * probe * obj[indy,:][:,indx] + \
                        (gamma/(2.0+gamma)) * phi[:,:,pos];
            old_phi_hat = scipy.fftpack.fft2(new_phi) / fnorm;
            phi_hat = magproj(old_phi_hat,amp_exp_norm[:,:,pos]);
            if not fmask is None:
                phi_hat = phi_hat * fmask[:,:,pos] + old_phi_hat * \
                            (fmask[:,:,pos] == 0).astype(fmask.dtype);
            phi[:,:,pos] = scipy.fftpack.ifft2(phi_hat) * fnorm;
        
        return Pty(phi,obj.copy(),probe.copy());





class PtychographyStats():
    """
    The Ptychography algorithms have special statistics. The necessary routines
    are collected in this class.
    """
    
    def __init__(self,config):
        """
        Initialization
        """
        self.Nx = config['Nx']; self.Ny = config['Ny'];
        self.dim = config['dim'];
        self.N_pie = config['N_pie'];
        self.norm_data = config['norm_data'];
        self.positions = config['positions'];
        self.object_support = config['object_support'];
        self.switch_object_support = config['switch_object_support_constraint'];
        self.probe = config['probe'];
        self.sample_plane = config['sample_plane'];
        self.fnorm = sqrt(self.Nx*self.Ny);
        self.amp_exp_norm = config['amp_exp_norm'];
        self.ptychography_prox = config['ptychography_prox'];
        self.fmask = config['fmask'];
        if not self.fmask is None:
            self.fmask = self.fmask.reshape((self.fmask.shape[0],
                                             self.fmask.shape[1],
                                             self.fmask.shape[2]*
                                             self.fmask.shape[3]));
    
    def objective(self,x1):
        """
        Get the value of the objective function (used in PALM)
        
        Parameters
        ----------
        x1 : Pty - Input
        
        Returns
        -------
        float - Objective function value of x1
        """
        obj_value = 0.0
        phi = x1.phi
        obj = x1.obj
        probe = x1.probe
        positions = self.positions
        
        rangeNy = numpy.arange(self.Ny,dtype=numpy.int)
        rangeNx = numpy.arange(self.Nx,dtype=numpy.int)
        
        for pos in range(self.N_pie):
            indy = rangeNy + positions[pos,0]
            indx = rangeNx + positions[pos,1]
            obj_value += scipy.linalg.norm(probe * obj[indy,:][:,indx] - \
                                                        phi[:,:,pos],'fro')**2
        return obj_value
    
    
    def change(self, x1, x2):
        """
        Computes the change for x1 and x2
        
        Parameters
        ----------
        x1 : Pty - Input 1
        x2 : Pty - Input 2
        
        Returns
        -------
        float - The change
        """
        norm_data = self.norm_data;
        
        changephi = 0.0;
        
        if self.Nx == 1 or self.Ny == 1:
            changephi = (scipy.linalg.norm(x1.phi-x2.phi,'fro')/norm_data)**2;
        else:
            for j in range(self.dim):
                changephi += (scipy.linalg.norm(x1.phi[:,:,j]-x2.phi[:,:,j], \
                                                            'fro')/norm_data)**2;
        if self.switch_object_support:
            changephi += (scipy.linalg.norm((x1.obj-x2.obj)*self.object_support,'fro')/norm_data)**2;
        else:
            changephi += (scipy.linalg.norm(x1.obj-x2.obj,'fro')/norm_data)**2;
        changephi += (scipy.linalg.norm(x1.probe-x2.probe,'fro')/norm_data)**2;
        return sqrt(changephi);
    
    
    def customerror(self, x1, x2, x3=None):
        """
        Computes the custom error for x1, x2 and x3
        
        Parameters
        ----------
        x1 : Pty - Input 1
        x2 : Pty - Input 2
        x3 : Pty, optional - Input 3
        
        Returns
        -------
        (rms_error_obj, rms_error_probe, change_phi, r_factor, objective) : tuple - The error
        """
        rangeNy = numpy.arange(self.Ny,dtype=numpy.int)
        rangeNx = numpy.arange(self.Nx,dtype=numpy.int)
        
        supported_plane = self.sample_plane * self.object_support
        
        rms_error_probe, diffphase, row_shift, col_shift = self._dftregistration( \
                scipy.fftpack.fft2(self.probe), scipy.fftpack.fft2(x2.probe))
        
        shifted_object = numpy.roll(numpy.roll(x2.obj,row_shift,axis=0), \
                                    col_shift,axis=1)
        supported_object = shifted_object * self.object_support
        gammaobj = numpy.sum(supported_plane * numpy.conj(supported_object)) / \
                    numpy.sum(shifted_object*numpy.conj(shifted_object)).real
        rms_error_obj = numpy.sum(numpy.abs(supported_plane - gammaobj * \
                        supported_object)**2) / numpy.sum(supported_plane * \
                        numpy.conj(supported_plane)).real
                        
        r_factor = 0.0
        for pos in range(self.N_pie):
            indy = rangeNy + self.positions[pos,0]
            indx = rangeNx + self.positions[pos,1]
            if self.fmask is None:
                r_factor += numpy.sum(numpy.abs(self.amp_exp_norm[:,:,pos] - \
                            numpy.abs(scipy.fftpack.fft2(x2.obj[indy,:][:,indx]* \
                            x2.probe)/self.fnorm)))
            else:
                r_factor += numpy.sum(numpy.abs(self.amp_exp_norm[:,:,pos] - \
                            numpy.abs(scipy.fftpack.fft2(x2.obj[indy,:][:,indx] * \
                            x2.probe)/self.fnorm))*self.fmask[:,:,pos])
        r_factor /= numpy.sum(self.amp_exp_norm)
        
        change_phi = 0.0;
        if self.Nx == 1 or self.Ny == 1:
            change_phi = (scipy.linalg.norm(x1.phi-x2.phi,'fro')/self.norm_data)**2
        else:
            for j in range(self.dim):
                change_phi += (scipy.linalg.norm(x1.phi[:,:,j]-x2.phi[:,:,j],'fro')/self.norm_data)**2
        
        if self.ptychography_prox == 'Thibault':
            x2 = x2.copy()
            x2.phi = x3.phi.copy()
        
        objective = self.objective(x2)
        
        return (rms_error_obj, rms_error_probe, change_phi, r_factor, objective)
        
    
    def _dftups(self, inp, nor=None, noc=None, usfac=1, roff=0, coff=0):
        """
        Upsampled DFT by matrix multiplications
        
        Parameters
        ----------
        inp : array_like - Input data
        nor, noc : int, optional - Number of pixels in the output unsampled DFT (Default = inp.shape)
        usfac : number, optional - Upsampling factor (Default = 1)
        roff, coff : int, optional - Row and column offsets (Default = 0)
            
        Returns
        -------
        ndarray - DFT
        """
        
        nr = inp.shape[0]; nc = inp.shape[1]
        if noc is None:
            noc = nc
        if nor is None:
            nor = nr
        kernc = numpy.exp((-2j*numpy.pi/(nc*usfac))*(scipy.fftpack.ifftshift( \
                        numpy.arange(nc))-(nc//2))*(numpy.arange(noc)-coff))
        kernr = numpy.exp((-2j*numpy.pi/(nr*usfac))*(scipy.fftpack.ifftshift( \
                        numpy.arange(nr))-(nr//2))*(numpy.arange(nor)-roff))
        return kernr*inp*kernc
        
    
    def _dftregistration(self, buf1ft, buf2ft, usfac=1):
        """
        Efficient subpixel image registration by crosscorrelation.
        
        Parameters
        ----------
        buf1ft : array_like - Fourier transform of the reference image
        buf2ft : array_like - Fourier transform of the image to register
        usfac : int, optional - Upsampling factor
        
        Returns
        -------
        error : number - Translation invariant normalized RMS error between f and g
        diffphase : number - Global phase difference between the two images
        row_shift : number - Pixel shifts between images (rows)
        col_shift : number - Pixel shifts between images (columns)
        
        References
        ----------
        .. [1] Manuel Guizar-Sicairos, Samuel T. Thurman, James R. Fienup,
           "Efficient subpixel image registration algorithms",
           Opt.Lett. 33, 156-158 (2008)
        """
        if usfac == 0:
            ccmax = numpy.sum(buf1ft * numpy.conj(buf2ft));
            rfzero = numpy.sum(buf1ft*numpy.conj(buf1ft));
            rgzero = numpy.sum(buf2ft*numpy.conj(buf2ft));
            error = 1.0 - ccmax * ccmax.conjugate() / (rfzero*rgzero);
            error = sqrt(abs(error));
            diffphase = atan2(ccmax.imag,ccmax.real);
            return error, diffphase, 0, 0;
        elif usfac == 1:
            m = buf1ft.shape[0]; n = buf1ft.shape[1];
            cc = scipy.fftpack.ifft2(buf1ft * numpy.conj(buf2ft));
            loc1 = numpy.argmax(numpy.abs(cc),axis=0);
            loc2 = numpy.argmax(numpy.abs(cc[loc1,numpy.arange(cc.shape[1])]));
            rloc = loc1[loc2];
            cloc = loc2;
            ccmax = cc[rloc,cloc];
            rfzero = numpy.sum(buf1ft*numpy.conj(buf1ft))/(m*n);
            rgzero = numpy.sum(buf2ft*numpy.conj(buf2ft))/(m*n);
            error = 1.0 - ((ccmax*ccmax.conjugate()) / (rgzero*rfzero));
            error = sqrt(abs(error));
            diffphase = atan2(ccmax.imag,ccmax.real);
            md2 = numpy.fix(m/2.0); nd2 = numpy.fix(n/2.0);
            row_shift = rloc - 1;
            col_shift = cloc - 1;
            if rloc > md2:
                row_shift -= m;
            if cloc > nd2:
                col_shift -= n;
            return error, diffphase, row_shift, col_shift;
        else:
            m = buf1ft.shape[0]; n = buf1ft.shape[1];
            mlarge = m*m; nlarge = n*n;
            cc = numpy.zeros((mlarge,nlarge));
            cc[m+1-numpy.fix(m/2.0):m+1+numpy.fix(m/2.0),n+1-numpy.fix(n/2.0):n+1+numpy.fix(n/2.0)] = \
                scipy.fftpack.fftshift(buf1ft) * numpy.conj(scipy.fftpack.fftshift(buf2ft));
            cc = scipy.fftpack.ifft2(scipy.fftpack.ifftshift(cc));
            loc1 = numpy.argmax(numpy.abs(cc),axis=0);
            loc2 = numpy.argmax(numpy.abs(cc[loc1,numpy.arange(cc.shape[1])]));
            rloc = loc1[loc2];
            cloc = loc2;
            ccmax = cc[rloc,cloc];
            m = cc.shape[0]; n = cc.shape[1];
            md2 = numpy.fix(m/2); nd2 = numpy.fix(n/2);
            row_shift = rloc - 1;
            col_shift = cloc - 1;
            if rloc > md2:
                row_shift -= m;
            if cloc > nd2:
                col_shift -= n;
            row_shift /= 2.0;
            col_shift /= 2.0;
            rg00 = None; rf00 = None;
            scale = 1.0/(md2*nd2*usfac*usfac);
            if usfac > 2:
                row_shift = round(row_shift*usfac)/usfac;
                col_shift = round(col_shift*usfac)/usfac;
                dftshift = numpy.fix(ceil(usfac*1.5)/2);
                cc = numpy.conj(self._dftups(buf2ft*numpy.conj( \
                        buf1ft),ceil(usfac*1.5),ceil(usfac*1.5), \
                        usfac,dftshift-row_shift*usfac, \
                        dftshift-col_shift*usfac)) * scale ;
                loc1 = numpy.argmax(numpy.abs(cc),axis=0);
                loc2 = numpy.argmax(numpy.abs(cc[loc1,numpy.arange(cc.shape[1])]));
                rloc = loc1[loc2];
                cloc = loc2;
                ccmax = cc[rloc,cloc];
                rg00 = self._dftups(buf1ft*numpy.conj(buf1ft),1,1,usfac)*scale;
                rf00 = self._dftups(buf2ft*numpy.conj(buf2ft),1,1,usfac)*scale;                
                rloc -= dftshift + 1;
                cloc -= dftshift + 1;
                row_shift += rloc/usfac;
                col_shift += cloc/usfac;
            else:
                rg00 = numpy.sum(numpy.real(buf1ft*numpy.conj(buf1ft)))/(m*n);
                rf00 = numpy.sum(numpy.real(buf2ft*numpy.conj(buf2ft)))/(m*n);
            error = 1.0 - ccmax*ccmax.conjugate() / (rg00*rf00);
            error = sqrt(abs(error));
            diffphase = atan2(ccmax.imag,ccmax.real);
            if md2 == 1:
                row_shift = 0;
            if nd2 == 1:
                col_shift = 0;
            return error, diffphase, row_shift, col_shift;


# We need these helper functions to circumvent some of the pitfalls of
# Python's multiprocessing
def _unwrap_do_block(arg, **kwarg):
    return Blocking_meta._do_block(*arg, **kwarg);

def _unwrap_do_view(arg, **kwarg):
    return Blocking_meta._do_view(*arg, **kwarg);


class Blocking_meta(Algorithm):
    """
    A meta-algorithm used if a blocking scheme is applied. Although this is
    probably not the most efficient implementation if blocking schemes, it
    can be applied independently of the algorithm.
    """
    
    def __init__(self, config):
        """
        Initialization
        """
        self.config = config.copy();
        
        self.subproblem_config = config.copy();
        self.subproblem_config['ignore_error'] = True;
        self.subsubproblem_config = self.subproblem_config.copy();
        self.subsubproblem_config['N_pie'] = 1;
        self.subsubproblem_config['dim'] = 1;
        
        self.within_blocks_scheme = config['within_blocks_scheme'];
        self.between_blocks_scheme = config['between_blocks_scheme'];
        
        self.ignore_error = config['ignore_error'];
        self.config['ptychography_prox'] = 'Blocking';
        self.custom_error = PtychographyStats(self.config);
        
        self.in_block = config['in_block'];
        self.num_blocks = numpy.amax(self.in_block)+1;
        
        self.block_iters = numpy.zeros((self.num_blocks,),dtype=numpy.int);
        
        self.algorithm = config['algorithm'];
        
        self.positions = config['positions'];
        self.amp_exp_norm = config['amp_exp_norm'];
        self.block_maxiter = config['block_maxiter'];
        self.block_tol = config['block_tol'];

        self.N_pie = config['N_pie'];
        
        self.pool = Pool();
    
    
    def _do_view(self, u, view, tol):
        """
        Creates and solves the subproblem for the current view
        
        Parameters
        ----------
        u : Pty - Input data
        view : array_like - Current view
        tol : number - Tolerance
        
        Returns
        -------
        sequence of array_like - A tuple of the form (phi,obj,probe)
        """
        cfg = self.subsubproblem_config.copy();
        
        cfg['positions'] = self.positions[[view],:];
        cfg['amp_exp_norm'] = self.amp_exp_norm[:,:,[view]];
        
        subsub_u = Pty(u.phi[:,:,[view]].copy(),u.obj.copy(),u.probe.copy());
        subsub_algorithm = self.algorithm(cfg);
        
        _u, result, _iters, _change, _gap = \
                                subsub_algorithm.run(subsub_u, tol, 1);
        
        return (result.phi, result.obj, result.probe);
        
    
    def _do_block(self, u, b, tol):
        """
        Creates and solves the subproblem for the current block
        
        Parameters
        ----------
        u : Pty - Input data
        b : int - Current block's index
        tol : number - Tolerance
        
        Returns
        -------
        sequence of array_like - A tuple of the form (phi,obj,probe)
        """
        in_block = self.in_block;
        block_maxiter = self.block_maxiter;
        within_blocks_scheme = self.within_blocks_scheme;
        N_pie = self.N_pie;
        positions = self.positions;
        amp_exp_norm = self.amp_exp_norm;
        
        B = (in_block == b);
        
        self.subproblem_config['positions'] = positions[B,:];
        self.subproblem_config['amp_exp_norm'] = amp_exp_norm[:,:,B];
        self.subproblem_config['N_pie'] = numpy.sum(B);
        self.subproblem_config['dim'] = self.subproblem_config['N_pie'];
        sub_u = Pty(u.phi[:,:,B].copy(),u.obj.copy(),u.probe.copy());
        
        sub_algorithm = self.algorithm(self.subproblem_config);
        
        sub_u_final = None;
        
        _block_iters = 0;
        
        if within_blocks_scheme == 'none':
            _u, sub_u_final, _iters, _change, _gap = \
                        sub_algorithm.run(sub_u, self.block_tol, block_maxiter);
            _block_iters = _iters;
        elif within_blocks_scheme == 'sequential':
            views_in_block = numpy.arange(N_pie)[B];
            
            sub_phi = numpy.empty((u.phi.shape[0], u.phi.shape[1], \
                            views_in_block.size),dtype=u.phi.dtype);
            
            for v in range(views_in_block.size):
                view = views_in_block[v];
                #self.subsubproblem_config['positions'] = \
                #                            self.positions[[view],:];
                #self.subsubproblem_config['amp_exp_norm'] = \
                #                        self.amp_exp_norm[:,:,[view]];
                #subsub_u = Pty(u.phi[:,:,[view]].copy(),u.obj.copy(),u.probe.copy());
                #
                #subsub_algorithm = self.algorithm(self.subsubproblem_config);
                #
                #_u, subsub_u_final, _iters, _change, _gap = \
                #                subsub_algorithm.run(subsub_u, tol, 1);
                subsub_u_final = self._do_view(u, view, tol);
                
                sub_phi[:,:,v] = subsub_u_final.phi[:,:,0];
            
            _block_iters = 1; # Sub-block iterations are limited to 1 
            
            sub_u_final = Pty(sub_phi, subsub_u_final.obj.copy(), \
                                          subsub_u_final.probe.copy());
        elif within_blocks_scheme == 'parallel':
            views_in_block = numpy.arange(N_pie)[B];
            views = range(views_in_block.size);
            
            arg = [(self, u, views_in_block[v], tol) for v in views];
            results = self.pool.map(_unwrap_do_view, arg);
            
            sub_phi = sum(t[0] for t in results);
            sub_obj = sum(t[1] for t in results) / views_in_block.size;
            sub_probe = sum(t[2] for t in results) / views_in_block.size;
            
            sub_u_final = Pty(sub_phi, sub_obj, sub_probe);            
            
        
        phi = numpy.zeros_like(u.phi);
        phi[:,:,B] = sub_u_final.phi;
        obj = sub_u_final.obj;
        probe = sub_u_final.probe;
        
        return (phi, obj, probe, _block_iters);
    
    
    def run(self, u, tol, maxiter):
        """
        Runs Blocking_meta
        """
        num_blocks = self.num_blocks;
        ignore_error = self.ignore_error;
        between_blocks_scheme = self.between_blocks_scheme;
        
        iters = 0;
        change = numpy.zeros(maxiter+1);
        change[0] = tol+1;
        
        if ignore_error == False:
            num_custom_errors = len(self.custom_error.customerror(u,u))
            custom_errors = numpy.zeros((num_custom_errors,maxiter+1));
        
        tmp_u = u.copy();
        
        while iters < maxiter and change[iters] > tol:
            iters += 1;
            
            if between_blocks_scheme == 'sequential':
                in_block = self.in_block;
                
                for b in range(num_blocks):
                    B = (in_block == b);
                    result = self._do_block(tmp_u, b, tol);
                    tmp_u.phi[:,:,B] = result[0][:,:,B];
                    tmp_u.obj = result[1];
                    tmp_u.probe = result[2];
                    self.block_iters[b] = result[3];
            elif between_blocks_scheme == 'parallel':    
                arg = [(self,tmp_u,b,tol) for b in range(num_blocks)];
                results = self.pool.map(_unwrap_do_block, arg);
                
                tmp_u.phi = sum([t[0] for t in results]);
                tmp_u.obj = sum([t[1] for t in results]) / num_blocks;
                tmp_u.probe = sum([t[2] for t in results]) / num_blocks;
                for b in range(num_blocks):
                    self.block_iters[b] = results[b][3];
            else:
                raise NotImplementedError('NYI')
            
            change[iters] = self.custom_error.change(u,tmp_u);
            
            if ignore_error == False:
                custom_errors[:,iters] = self.custom_error.customerror(u,tmp_u);
            
            u = tmp_u.copy();
        
        change = change[1:iters+1];
        custom_errors = numpy.empty((5,0)) if ignore_error else custom_errors[:,1:iters+1];
        
        return u, u.copy(), iters, change, custom_errors;
    
    
    # More multiprocessing workarounds
    def __getstate__(self):
        """
        a multiprocessing workaround
        """
        self_dict = self.__dict__.copy();
        del self_dict['pool'];
        return self_dict;

    def __setstate__(self, state):
        """
        a multiprocessing workaround
        """
        self.__dict__.update(state);



class PALM(Algorithm):
    """
    Proximal alternating (linearized) minimization algorithm
    
    """
    
    def __init__(self, config):
        """
        Parameters
        ----------
        config : dict        
            Dictionary containing the problem configuration.
            It must contain the following mappings:
            
            projector_sequence: sequence of Proxoperator
                 Sequence of prox operators to be used. (The classes, no instances)
            Nx: int
                Row-dim of the product space elements
            Ny: int
                Column-dim of the product space elements
            Nz: int
                Depth-dim of the product space elements
            dim: int
                Size of the product space
            norm_data: number
                ?
            ignore_error: boolean
                Whether to ignore errors
        """
        self.projs = [p(config) for p in config['projector_sequence']];
        
        self.ignore_error = config['ignore_error'];
        self.custom_error = PtychographyStats(config);
    
    
    def run(self, u, tol, maxiter):
        """
        Runs the algorithm one some input data
        
        Parameters
        ----------
        u : array_like - Input data
        tol : number - Tolerance
        maxiter : int - Maximum number of iterations
        
        Returns
        -------
        u : array_like - Result
        u_final : array_like - Result
        iters : int - Number of iterations the algorithm performed
        change : array_like - The percentage change in the norm
        gap : array_like - Squared gap distance normalized by the magnitude constraint
        """
        projs = self.projs;
        ignore_error = self.ignore_error;
        
        iters = 0;
        change = numpy.zeros(maxiter+1);
        change[0] = tol+1;
        
        if ignore_error == False:
            num_custom_errors = len(self.custom_error.customerror(u,u))
            custom_errors = numpy.zeros((num_custom_errors,maxiter+1));
        
        tmp_u = u;
        
        while iters < maxiter and change[iters] > tol:
            iters += 1;
            
            for p in projs:
                tmp_u = p.work(tmp_u);
            
            change[iters] = self.custom_error.change(u,tmp_u);
            
            if self.ignore_error == False:
                custom_errors[:,iters] = self.custom_error.customerror(u,tmp_u);
                
            u = tmp_u;
        
        change = change[1:iters+1];
        custom_errors = numpy.empty((5,0)) if ignore_error else custom_errors[:,1:iters+1];
        
        return tmp_u, tmp_u.copy(), iters, change, custom_errors;


class RAAR_PALM(Algorithm):
    """
    RAAR implementation that produces error statistics comparable to PALM
    
    """
    
    def __init__(self,config):
        """
        Parameters
        ----------
        config : dict        
            Dictionary containing the problem configuration.
            It must contain the following mappings:
        
            projector_sequence: sequence of ProxOperator
                Sequence of prox operators to be used. (The classes, no instances).
                For RAAR_PALM it must contain 2 ProxOperators.
            Nx: int
                Row-dim of the product space elements
            Ny: int
                Column-dim of the product space elements
            Nz: int
                Depth-dim of the product space elements
            dim: int
                Size of the product space
            beta0: number
                Starting relaxation parmater
            beta_max: number
                Maximum relaxation parameter
            beta_switch: int
                Iteration at which beta moves from beta0 -> beta_max
            norm_data: number
                ?
            ignore_error: boolean
                Whether to ignore errors
        """
        self.beta0 = config['beta0'];
        self.beta_max = config['beta_max'];
        self.beta_switch = config['beta_switch'];
        
        self.projs = [p(config) for p in config['projector_sequence']];
        
        self.ignore_error = config['ignore_error'];
        self.custom_error = PtychographyStats(config);
    
    
    def run(self, u, tol, maxiter):
        """
        Runs the algorithm one some input data
        
        Parameters
        ----------
        u : array_like - Input data
        tol : number - Tolerance
        maxiter : int - Maximum number of iterations
        
        Returns
        -------
        tmp_u1 : array_like - Result
        tmp_u1.copy() : array_like - Copy of result
        iters : int - Number of iterations the algorithm performed
        change : array_like - The percentage change in the norm
        custom_errors : array_like - Squared gap distance normalized by the magnitude constraint
        """
        projs = self.projs;
        ignore_error = self.ignore_error;
        
        iters = 0
        change = numpy.zeros(maxiter+1)
        change[0] = tol+1
        
        if ignore_error == False:
            num_custom_errors = len(self.custom_error.customerror(u,u,u))
            custom_errors = numpy.zeros((num_custom_errors,maxiter+1))
        
        tmp_u1 = projs[0].work(u)
        
        while iters < maxiter and change[iters] > tol:
            iters += 1
            tmp = exp((-iters/self.beta_switch)**3)
            beta = (tmp*self.beta0) + ((1-tmp)*self.beta_max)
            
            tmp_u3 = tmp_u1.copy()
            tmp_u3.phi = 2*tmp_u1.phi - u.phi
            
            tmp_u4 = projs[1].work(tmp_u3)
            tmp_u = tmp_u4.copy()
            tmp_u.phi = (beta*(2*tmp_u4.phi-tmp_u3.phi)+(1-beta)*tmp_u3.phi+u.phi)/2
            
            tmp_u2 = projs[0].work(tmp_u)
            
            change[iters] = self.custom_error.change(tmp_u1,tmp_u2)
            
            if self.ignore_error == False:
                custom_errors[:,iters] = self.custom_error.customerror(tmp_u1,tmp_u2,tmp_u4)
            
            u = tmp_u
            tmp_u1 = tmp_u2
        
        change = change[1:iters+1]
        custom_errors = numpy.empty((5,0)) if ignore_error else custom_errors[:,1:iters+1]
        
        return tmp_u1, tmp_u1.copy(), iters, change, custom_errors




class Ptychography(Problem):
    """
    Ptychography Problem
    
    References
    ----------
    .. [1] R. Hesse, D.R. Luke, S. Sabach, and M.K. Tam, "Proximal
       Heterogeneous Block Implicit-Explicit Method and Application to Blind
       Ptychographic Diffraction Imaging", SIAM J. Imaging Sciences,
       8(1):426–457, 2015.
    
    """
    config = {
        'sim_data_type':'gaenseliesel',
        'Nx':64,
        'Ny':64,
        'Nz':1,
        'scan_stepsize':3.5e-7,
        'nx':25,
        'ny':25,
        'sample_area_center':(0.5,0.5),
        'noise':None,
        'switch_probemask':True,
        'probe_mask_gamma':1,
        'rmsfraction':0.5,
        'scan_type':'raster',
        'switch_object_support_constraint':True,
        'probe_guess_type':'circle',
        'object_guess_type':'random',
        'warmup':True,
        'warmup_alg':'PALM',
        'warmup_maxiter':2,
        'algorithm':'PALM',
        'maxiter':2, # This is usually 300. I reduced it to be able to test it.
        'tol':0.625,
        'ignore_error':False,
        'ptychography_prox':'Rodenburg',
        'blocking_switch':True,
        'blocking_scheme':'divide',
        'between_blocks_scheme':'sequential',
        'within_blocks_scheme':'parallel',
        'block_rows':2,
        'block_cols':2,
        'block_maxiter':1,
        'block_tol':1.0,
        'rodenburg_inner_it':1,
        'beta0':1,
        'beta_max':1,
        'beta_switch':30,
        'bs_mask':None,
        'bs_factor':1,
        'fmask':None,
        'overrelax':1
    }
    
    
    def __init__(self,new_config={}):
        """
        Initialization of Ptychography problem
        
        Parameters
        ----------
        config : dict, optional - Dictionary containing the problem configuration. If unspecified, Ptychography.config is used.
        """
        self.config.update(new_config)
        
        #self.config['algorithm'] = getattr(Algorithms, self.config['algorithm'])
        self.config['algorithm'] = globals()[self.config['algorithm']]
        self.config['warmup_alg'] = globals()[self.config['warmup_alg']]
        
        if self.config['ptychography_prox'] == 'PALM':
            self.config['projector_sequence'] = [P_PHeBIE_probe,
                                                 P_PHeBIE_object,
                                                 P_PHeBIE_phi];
            self.config['warmup_projsequence'] = [P_PHeBIE_object,
                                                  P_PHeBIE_phi];
        elif self.config['ptychography_prox'] == 'PALMregPhiPtwise':
            self.config['projector_sequence'] = [P_PHeBIE_probe_ptwise,
                                                 P_PHeBIE_object_ptwise,
                                                 P_PHeBIE_phi_regularized];
            self.config['warmup_projsequence'] = [P_PHeBIE_object_ptwise,
                                                  P_PHeBIE_phi_regularized];
        elif self.config['ptychography_prox'] == 'Rodenburg':
            self.config['projector_sequence'] = [P_rodenburg];
            self.config['warmup_projsequence'] = [P_rodenburg_probe_fixed];
        elif self.config['ptychography_prox'] == 'Thibault':
            self.config['projector_sequence'] = [P_thibault_op,P_thibault_f];
            self.config['warmup_projsequence'] = [P_thibault_o,P_thibault_f];
            self.config['algorithm'] = RAAR_PALM;
            self.config['warmup_alg'] = RAAR_PALM;
#            print ('Operator set \'Thibault\' requires the use of RAAR_PALM')
        elif self.config['ptychography_prox'] == 'Thibault_AP':
            self.config['projector_sequence'] = [P_thibault_op,P_thibault_f];
            self.config['warmup_projsequence'] = [P_thibault_o,P_thibault_f];
        else:
            raise NotImplementedError ( \
                'Operator set \'%s\' not (yet) supported'%self.config['ptychography_prox'])
        
        self.config['N_pie'] = self.config['nx']*self.config['ny']
        self.config['dim'] = self.config['N_pie']; # warum zwei Parameter mit gleichem Wert vorhalten?
        self.config['norm_data'] = 1
    

    def _gaussu(self, rho, z, l, wnull):
        """
        Calculates the electric field component of a Gaussian beam
        
        Parameters
        ----------
        rho : array_like - Input image
        z : int - Position on the Z axis
        l : int/float - Wave length
        wnull : number - Waist radius
        
        Returns
        -------
        array_like - Field amplitude
        """
        k = 2.0*numpy.pi/l;
        anull = 1.0;
        znull = (wnull**2) * numpy.pi/l;
        w = wnull * sqrt(1 + (z**2)/(znull**2));
        
        if z == 0:
            r = numpy.Inf;
        else:
            r = z * (1 + (znull**2)/(z**2));
        
        u = anull * (1/r + 2j/(k*w*w)) * numpy.exp(1j*k*z) *  \
            numpy.exp((1j*k/(2*r) - 1/(w*w)) * (rho**2));
        return u;
    
    
    def _prop_nf_pa(self,im,l,delta_z,dx,dy):
        """
        Propagates a 2D-wave field within a paraxial approximation into the
        optical near field.
        
        Parameters
        ----------
        im : array_like - Incident 2D wave field
        l : number - Wave length
        delta_z : number - Propagation distance
        dx, dy : int - Pixel widths in horizontal and vertical direction
        
        Returns
        -------
        array_like - Propagated wave field
        
        """
        Nx = im.shape[1]; Ny = im.shape[0];
        k = 2*numpy.pi/l;
        
        if delta_z == 0:
            alpha_th = numpy.Inf;
        else:
            alpha_th = (4/numpy.pi * l/delta_z)**0.25;
        alpha = max([atan(Nx/dx*(2*delta_z)),atan(Ny/dy*(2*delta_z))]);
        if alpha > alpha_th:
                raise Exception('Small-angle approximation is violated')
        
        px = abs(l*delta_z/(Nx*dx*dx));
        py = abs(l*delta_z/(Ny*dy*dy));
        f = max([px,py]);
        
        beta = 1.0;
        
        if f > 1:
            # This is not yet tested
            print ('Probe-FOV is enlarged by factor %3.2f for propagation' % beta*f)
            Nx_l = ceil(beta*f)*Nx;
            dqx = 2*numpy.pi/(Nx_l*dx);
            Ny_l = ceil(beta*f)*Ny;
            dqy = 2*numpy.pi/(Ny_l*dy);
            
            rangeNx = numpy.arange(Nx_l);
            rangeNy = numpy.arange(Ny_l);
            
            alpha_th = (4/numpy.pi*l/delta_z)**0.25;
            alpha = max([atan(Nx/dx*(2*delta_z)),atan(Ny/dy*(2*delta_z))]);
            if alpha > alpha_th:
                raise Exception('Small-angle approximation is violated')
            
            Qx,Qy = numpy.meshgrid((rangeNx-(Nx_l//2))*dqx,(rangeNy-(Ny_l//2))*dqy);
            kappa = -0.5*scipy.fftpack.ifftshift(Qx**2 + Qy**2)/k;
            
            Cx = (Nx_l//2)+1; Cy = (Ny_l//2)+1;
            im_l = numpy.zeros(Ny_l,Nx_l);
            
            for i in range(Cy-(Ny//2)-1):
                im_l[i,Cx-(Nx//2)-1:Cx+ceil(Nx/2.0)-1] = im[0,:];
            for i in range(Cy+ceil(Ny/2.0)-1,Ny_l):
                im_l[i,Cx-(Nx//2)-1:Cx+ceil(Nx/2.0)-1] = im[-1,:];
            for i in range(Cx-(Nx//2)-1):
                im_l[Cy-(Ny//2)-1:Cy+ceil(Ny/2.0)-1,i] = im[:,0];
            for i in range(Cx+ceil(Nx/2.0)-1,Nx_l):
                im_l[Cy-(Ny//2)-1:Cy+ceil(Ny/2.0)-1,i] = im[:,-1];
            im_l[0:Cy-(Ny//2)-1,0:Cx-(Nx//2)-1] = im[0,0];
            im_l[0:Cy-(Ny//2)-1,Cx+ceil(Nx/2.0)-1:] = im[0,-1];
            im_l[Cy+ceil(Ny/2.0)-1:Ny_l,0:Cx-(Nx//2)-1] = im[-1,0];
            im_l[Cy+ceil(Ny/2.0)-1:Ny_l,Cx+ceil(Nx/2.0)-1:] = im[-1,-1];
            
            im_l[Cy-(Ny//2)-1:Cy+ceil(Ny/2.0)-1,Cy-(Ny//2)-1:Cy+ceil(Ny/2.0)-1] = im;
            
            Im_l = scipy.fftpack.fftshift(scipy.fftpack.ifft2( \
                    scipy.fftpack.fft2(im_l) * numpy.exp(1j*kappa*delta_z)));
            return Im_l[Cy-(Ny//2)-1:Cy+ceil(Ny/2.0)-1,Cy-(Ny//2)-1:Cy+ceil(Ny/2.0)-1];
        else:
            dqx = 2*numpy.pi/(Nx*dx); dqy = 2*numpy.pi/(Ny*dy);
            rangeNx = numpy.arange(Nx);
            rangeNy = numpy.arange(Ny);
            Qx,Qy = numpy.meshgrid((rangeNx-(Nx//2))*dqx,(rangeNy-(Ny//2)*dqy));
            kappa = -0.5*scipy.fftpack.ifftshift(Qx**2 + Qy**2)/k;
            return scipy.fftpack.fftshift(scipy.fftpack.ifft2(scipy.fftpack.fft2( \
                scipy.fftpack.ifftshift(im)) * numpy.exp(1j*kappa*delta_z)));
    
    
    def _im2hsv(self,im,th):
        """
        method
        """
        abs_im = numpy.abs(im); abs_im -= numpy.amin(abs_im);
        abs_im /= numpy.amax(abs_im+numpy.spacing(1.));
        
        imhsv = numpy.zeros((im.shape[0],im.shape[1],3));
        
        imhsv[:,:,0] = numpy.remainder(numpy.angle(im),2*numpy.pi) / (2*numpy.pi);
        imhsv[:,:,1] = 1;
        tmp = (1.0 + numpy.sign(abs_im-th)) / 2.0;
        tmp = numpy.clip((abs_im - tmp*abs_im + tmp*th) / th, 0, 1);
        imhsv[:,:,2] = tmp;
        
        return imhsv;
    
    
    def _get_ptychography_data(self):
        """
        Generates an experiment on-the-fly. If you want to load existing
        data or make your own experiments, you can overload this.
        
        Returns
        -------
        d1x : number - Pixel dimension
        d1y : number - Pixel dimension
        positions : array_like
        probe : array_like
        sample_plane : array_like
        i_exp : array_like
                 
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
        w0 = fwhm/(2*sqrt(log(2)));
        
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
                assert( float(int(origin_x0-1)) == (origin_x0-1))
                assert( float(int(origin_y0-1)) == (origin_y0-1))
                i_exp[:,:,akk,bkk] = numpy.abs(numpy.rot90(numpy.roll( \
                    numpy.roll(scipy.fftpack.fft2(esw),int(origin_x0-1),axis=0), \
                    int(origin_y0-1),axis=1),-2))**2;
                
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
        nx = self.config['nx'];
        ny = self.config['ny'];
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
            in_block = numpy.arange(N_pie) % cols;
            for i in range(nx):
                in_block[i*nx:(i+1)*nx] += (i % rows) * cols;
        elif blocking_scheme == 'split':
            in_block = numpy.arange(N_pie) // (N_pie//(rows*cols));
        else:
            raise NotImplementedError('NYI')
        self.config['in_block'] = in_block;
    
    
    def _presolve(self):
        """
        Prepares argument for actual solving routine
        """
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
                    if sqrt((i-x_c)**2 + (j-y_c)**2) < radius:
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
                                sqrt(numpy.size(probe_guess));
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
            object_guess = numpy.ones_like(sample_plane);
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
        """
        Runs the warmup-algorithm if configured and
        starts the actual solving process
        to solve the given ptychography problem
        """
        if self.config['warmup'] == True:
#            print ('Starting warm-up algorithm')
            u, self.u_warm, iters, gap, change = \
                self.warmup_alg.run(self.u,self.config['tol'],self.config['warmup_maxiter']);
#            print ('End warm-up algorithm')
        else:
#            print ('Skipping warm-up algorithm')
            self.u_warm = self.u;
        
#        print ('Start algorithm')
        u, self.u_final, self.iters, self.change, self.custom_errors = \
                self.algorithm.run(self.u_warm,self.config['tol'],self.config['maxiter']);
#        print ('End algorithm')
    
    
    def _postsolve(self):
        """
        Processes the solution and generates the output
        """
        # probably nothing to do except for displaying things
        # hence, see method 'show'
        return
        
    def show(self):
        """
        This procedure generates several figures visualizing the result.
        """
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
        pyplot.imshow(colors.hsv_to_rgb(self._im2hsv(probe,1.0)));
        ax = pyplot.subplot(1,2,2);
        ax.title.set_text('Exact probe');
        pyplot.imshow(colors.hsv_to_rgb(self._im2hsv(self.config['probe'],1.0)));
        
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
    


class Ptychography_NTT_01_26210(Ptychography):
    """
    Ptychography NTT 01 26210
    """

    default_config = {
        'Nx':192,
        'Ny':192,
        'Nz':1,
        'scan_stepsize':3.5e-7,
        'nx':26,
        'ny':26,
        'sample_area_center':(0.5,0.5),
        'noise':None,
        'switch_probemask':False,
        'probe_mask_gamma':1,
        'rmsfraction':0.5,
        'scan_type':'raster',
        'switch_object_support_constraint':False,
        'probe_guess_type':'robin_initial',
        'object_guess_type':'constant',
        'warmup':True,
        'warmup_alg':PALM,
        'warmup_maxiter':0,
        'algorithm':PALM,
        'maxiter':2, # This is usually 500. I reduced it to be able to test it.
        'tol':-1,
        'ignore_error':False,
        'ptychography_prox':'Rodenburg',
        'blocking_switch':False,
        'blocking_scheme':'single_view',
        'between_blocks_scheme':'averaged',
        'within_blocks_scheme':'none',
        'block_rows':2,
        'block_cols':2,
        'block_maxiter':1,
        'block_tol':1.0,
        'rodenburg_inner_it':1,
        'beta0':1,
        'beta_max':1,
        'beta_switch':30,
        'bs_mask':None,
        'bs_factor':1,
        'overrelax':1
    }
    
    def __init__(self,config=default_config):
        """
        Initialization based on :class:`Ptychography`
        """
        Ptychography.__init__(self,config);
    
    def _get_ptychography_data(self):
        """
        Overloads the corresponding method of :class:`Ptychography`
        to create other experiment data
        """
        data = scipy.io.loadmat('../inputdata/data_NTT_01_26210_192x192.mat');
        d1y = data['d1y'][0][0];
        d1x = data['d1x'][0][0];
        i_exp = data['I_exp'];
        l = data['lambda'][0][0];
        z01 = data['z01'][0][0];
        self.config['fmask'] = data['fmask'].astype(numpy.float);
        data = scipy.io.loadmat('../inputdata/reconstruction_data_NTT_01_26210_192x192.mat');
        probe = data['Probe'];
        sample_plane = data['object'];
        
        data = scipy.io.loadmat('../inputdata/positions_NTT_01_26210.mat')
        positions = data['positions'].astype(numpy.float);
        positions[:,0] /= d1y;
        positions[:,0] -= numpy.amin(positions[:,0]);
        positions[:,1] /= d1x;
        positions[:,1] -= numpy.amin(positions[:,1]);
        positions[numpy.modf(positions)[0] == 0.5] += 0.01; # numpy.round is weird
        positions = numpy.round(positions).astype(numpy.int);
        
        self.config['trans_max_true'] = 1.0;
        self.config['trans_min_true'] = 0.0;
        
        return d1x, d1y, l, z01, positions, probe, sample_plane, i_exp;

