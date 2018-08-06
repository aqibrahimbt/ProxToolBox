from pathlib import Path
import numpy as np
from numpy import exp, sqrt, log2,log10, floor, unravel_index, argmax, zeros
from scipy.io import loadmat
from numpy.fft import fftshift, ifft2, fft2
from proxtoolbox.Utilities.mypoissonrnd import mypoissonrnd
from proxtoolbox.Utilities.gaussian import gaussian
from matplotlib.pyplot import subplots, show

def Goettingen_data_processor(config):
#  Parameters of the forward problem


    my_file = Path("../InputData/Phase/cell256.mat")
    if not(my_file.is_file()):
        print("Phase input data is missing.") 
        if GetData.query_yes_no("Do you want to download the phase input data?"):
            urllib.request.urlretrieve("http://num.math.uni-goettingen.de/~r.luke/tmp/Phase.tar.gz","../InputData/Phase.tar.gz", reporthook=GetData.dlProgress)
            print("\nExtracting data...")
            tar = tarfile.open("../InputData/Phase.tar.gz", "r:gz")
            tar.extractall("../InputData/")
            tar.close()
    if not(my_file.is_file()):
            print('***************************************************************************************')
            print('* Input data still missing.  Please try automatic download again or manually download *') 
            print('*    http://num.math.uni-goettingen.de/data/Phase.tar.gz                              *')
            print('* Save and unpack the Phase.tar.gz datafile in the                                    *')
            print('*    ProxMatlab/InputData subdirectory                                                *')
            print('***************************************************************************************')


    experiment = config['experiment']
    
    if(experiment =='CDI'):
            
        print('Loading data file CDI_intensity')
        f = loadmat('../InputData/Phase/CDI_intensity.mat')
        # diffraction pattern
        dp = f['intensity']['img'][0,0]
        orig_res = dp.size # actual data size
        step_up = np.ceil(log2(config['Nx'])-log2(orig_res))
        workres = 2**(step_up)*2**(floor(log2(orig_res)))#desired array size
    
        N=int(workres)
    
        ## center data and use region of interest around center
        #central pixel
        #find max data value
        argmx = unravel_index(argmax(dp), dp.shape)
    
        Xi = argmx[0]+1
        Xj = argmx[1]+1
    
        #get desired roi:
        #necessary conversion
        Di = N/2-(Xi-1)
        Dj = N/2-(Xj-1)
    
        # roi around central pixel
        Ni = 2*(Xi+Di*(Di<0)-1)
        Nj = 2*(Xj+Dj*(Dj<0)-1)
    
        tmp = zeros((N,N))
    
        tmp[int(Di*(Di>0)):int(Di*(Di>0)+Ni),int(Dj*(Dj>0)):int(Dj*(Dj>0)+Nj)] = dp[int(Xi-Ni/2)-1:int(Xi+Ni/2-1),int(Xj-Nj/2)-1:int(Xj+Nj/2-1)]
    
        M=(fftshift(tmp))**(.5)
        # M=tmp.^(.5)
        ## define support
        DX = np.ceil(N/7)
        S = zeros((N,N))
        S[int(N/4+1+DX)-1:int(3/4*N-1-DX),int(N/4+1+DX)-1:int(3/4*N-1-DX)] = 1
    
        config['magn'] = 1 # magnification factor
        if config['fresnel_nr']*config['magn'] <=  2*np.pi*sqrt(config['Nx']*config['Ny']):
            config['use_farfield_formula'] = 1
            print('Using farfield formula.')
        else:
            config['use_farfield_formula'] = 0
            print('Using nearfield formula.')
    
    
        # config['data_ball']=config['Nx']*config['Ny']*data_ball
        config['rt_data']=M
        # standard for the main program is that 
        # the data field is the magnitude SQUARED
        # in Luke.m this is changed to the magnitude.
        config['data']=M**2
        config['norm_rt_data']=np.linalg.norm(ifft2(M),'fro')
        config['data_zeros'] = np.where(M==0)
        config['support_idx'] = np.nonzero(S)
        config['product_space_dimension'] = 1
    
        # use the abs_illumination field to represent the 
        # support constraint.
        config['abs_illumination'] = S 
        config['supp_phase'] = []
    
        # initial guess
        config['u_0'] = S*np.random.rand(N,N)
        config['u_0'] = config['u_0']/np.linalg.norm(config['u_0'],'fro')*config['norm_rt_data'] 
        
        
        f, (ax1, ax2) = subplots(1, 2)
    
        im=ax1.imshow(log10(dp))
        f.colorbar(im, ax=ax1)
        ax1.set_title('Far field data')
        im=ax2.imshow(config['abs_illumination'])
        ax2.set_title('Support constraint')
        show()

    elif (experiment == 'Krueger'):
        # Sample: NTT-AT, ATN/XREESCO-50HC
        # 500 nm thick Ta structure: amplitude transmission 0.9644,
        # phase shift 0.4rad at E = 17.5 keV
        # parameters see below

        # empty waveguide (WG) beam
        WG = loadmat('../InputData/Phase/WG_beam.mat')
        WG = WG['WG']

        # hologram
        print('Loading data hologram_not-normalized.mat')
        I_exp = loadmat('../InputData/Phase/hologram_not-normalized.mat')
        I_exp = I_exp['I_exp']

        I_exp[np.isnan(I_exp)] = 1


        # The following is NOT used because it is not clear how the ``normalized" hologram 
        # is obtained.  I suspect it is obtained by simply dividing out the empty beam
        # data in the obervation plane.  But this is incorrect.  Despite the 
        # efforts of the Hohage team to get around this error by ad hoc regularization, 
        # we take the following exact approach:  back propagate the empty beam, and 
        # divide this from the unknown object in the object plane.  This approach is true 
        # to the imaging model and does not introduce any approximations. It does, however, 
        # expose the fact that we do not know the phase of the beam.   
        # clear I_exp
        # load 'data/hologram.mat'
        ## single image reconstruction


        # total number of elements
        N = I_exp.size

        #number of pixels
        Ny = I_exp.shape[0]; Nx = I_exp.shape[1]
        config['Ny']=Ny; config['Nx']=Nx

        ##########################
        # Experimental parameters
        ##########################

        # energy in keV
        E = 17.5

        # wavelength [m]
        lambd = 12.398/E*1E-10
        k = 2*np.pi/lambd

        # distance source-sample [m]
        z1 = 7.48e-3

        # distance sample-detector [m]
        z2 = 3.09

        # effective distance of detector
        z_eff = z1*z2/(z1+z2)

        # effective magnification factor
        M = (z1+z2)/z1

        #pixel size in detector plane
        dx_det = 55e-6
        dy_det = 55e-6

        # magnified coordinate system in detector plane
        # [X_det,Y_det] = meshgrid(dx_det*[0:1:Nx-1],dy_det*[0:1:Ny-1]);

        # demagnifiedpixel size in sample plane
        dx = dx_det/M
        dy = dy_det/M

        # coordinate system in sample plane
        # [X,Y] = meshgrid(dx*((1:1:Nx)-floor(Nx/2)-1),dy*((1:1:Ny)-floor(Ny/2)-1));

        # magnified coordinate system in detector plane
        # [X_det,Y_det] = meshgrid(dx_det*((1:1:Nx)-floor(Nx/2)-1),dy_det*((1:1:Ny)-floor(Ny/2)-1));

        # grid conversion in q-space
        dqx = 2*np.pi/(Nx*dx)
        dqy = 2*np.pi/(Ny*dy)

        [Qx,Qy] = np.meshgrid(dqx*(range(1,Nx+1)-np.floor(Nx/2)-1),dqy*(range(1,Ny+1)-np.floor(Ny/2)-1))

        ##################################
        # Prepare data from ProxToolbox:
        ##################################
        # Fresnel propagator:
        config['use_farfield_formula']=False
        config['Nx']=Nx
        config['Ny']=Ny
        config['fresnel_nr'] = 1*2*np.pi*Nx # Fresnel number
        config['magn']=1 # magnification factor (should this be M above, or 
        # is it okay since the demagnified pixel size is used?)

        kappa = np.sqrt(k**2-(Qx**2+Qy**2))
        config['FT_conv_kernel'] = fftshift(exp(-1j*kappa*z_eff))
        config['beam'] = abs(ifft2(fft2(sqrt(WG*7.210116465530644e-04))/config['FT_conv_kernel']))
        # the scaling factor of the waveguide array above was 
        # reverse engineered from the scaling that was apparently
        # used in the preparation of the data in the file hologram.mat
        # I don't want to use this hologram as the data, preferring instead
        # to keep the experimental data as close to the actual observation
        # as possible.  Also, I am disagreeing with what my colleagues in 
        # physics think is the correct way to compensate for a structured 
        # empty beam.  According to my reverse engineering, it appears that 
        # they have divided the EMPTY BEAM MEASUREMENT WG from the object
        # IN THE OBJECT PLANE.  It is true that you need to do the empty beam 
        # correction in the object plane, though not with the empty beam, but 
        # rather the reverse propagated empty beam, as I have done above. The
        # reason being that the empty beam is measured in the measurement plane, 
        # so it does not make sense to divide the object&beam in the object 
        # plane by the empty beam in the measurement plane -  you should divide
        # by the empty beam propagated back to the object plane.  The difficulty in 
        # this is that we do not know the phase of the empty beam.  The above 
        # formulation sets the phase of the beam to zero at the object plane. 
        # The beam is saved as the square root of the empty beam to conform with 
        # the projeciton operations working on the amplitude instead of the
        # magnitude.			

        # problem data
        config['data'] = I_exp
        config['rt_data'] = sqrt(I_exp)
        config['data_zeros'] = config['data']==0
        config['norm_rt_data'] =sqrt(sum(sum(config['data'])))
        # use the abs_illumination field to represent the 
        # support constraint.
        config['abs_illumination'] = np.ones(I_exp.shape)
        config['support_idx'] = config['abs_illumination'].nonzero()
        config['product_space_dimension'] = 1
        config['supp_phase'] = []

        # start values
        config['u_0']=ifft2(fft2(config['rt_data']*config['magn'])/config['FT_conv_kernel'])/config['beam']

        '''
        figure(1);
        subplot(2,2,1)
        imagesc(Nx*dx*linspace(0,1,Ny)*1E6,Nx*dx*linspace(0,1,Nx)*1E6,config['beam);
        axis equal;
        axis tight;
        xlabel('x [\mum]');
        ylabel('y [\mum]');
        colormap(gray); colorbar; title('empty beam - detector plane')

        subplot(2,2,2)
        imagesc(Nx*dx*linspace(0,1,Ny)*1E6,Nx*dx*linspace(0,1,Nx)*1E6,config['data);
        axis equal;
        axis tight;
        xlabel('x [\mum]');
        ylabel('y [\mum]');
        colormap(gray); colorbar; title('uncorrected near field observation')


        subplot(2,2,3)
        imagesc(Nx*dx*linspace(0,1,Ny)*1E6,Nx*dx*linspace(0,1,Nx)*1E6,(abs(config['u_0)));
        axis image; colormap(gray); colorbar;
        xlabel('x [\mum]');
        ylabel('y [\mum]');
        title('initial guess amplitude');

        # # pattern
        subplot(2,2,4);
        imagesc(Nx*dx*linspace(0,1,Ny)*1E6,Nx*dx*linspace(0,1,Nx)*1E6,(angle(config['u_0))');
        axis image; colormap(gray); colorbar;
        xlabel('x [\mum]');
        ylabel('y [\mum]');
        title('initial guess phase');
        caxis([-0.9 -0.4]);
        '''
        
    elif(experiment == 'dictyM103_stx6_600frames'):
        # project stx6 
        # folder_base: 'Z:\analysis\petraIII\run7\framework\dicty_M103'
        # DRL 07.06.2017:  I don't think the parameter values are correct. 
        # hologram
        print('Loading data dictyM103_stx6_600frames')
        data = loadmat('../InputData/Phase/dictyM103_stx6_600frames.mat')
        P = data['P'][0,0]
        
        # total number of elements
        N = data['DP'].size
    
        #number of pixels
        Ny = data['DP'].shape[0]; Nx = data['DP'].shape[1]
        config['Ny']=Ny; config['Nx']=Nx

        ##########################
        # Experimental parameters
        ##########################
    
        # energy in keV
        E = P['exp']['E'].astype(np.float32)
    
        # wavelength [m]
        lambd = P['exp']['lambda'].astype(np.float32)
        k = P['exp']['k'].astype(np.float32)
    
        # effective distance of detector
        z_eff = P['exp']['z_eff'].astype(np.float32)
    
        # effective magnification factor
        M = P['exp']['M'].astype(np.float32)
    
        #pixel size in detector plane
        dx_det = P['det']['dx'].astype(np.float32)
        dy_det = P['det']['dx'].astype(np.float32)
    
        # magnified coordinate system in detector plane
        # [X_det,Y_det] = meshgrid(dx_det*[0:1:Nx-1],dy_det*[0:1:Ny-1]);
    
        # demagnifiedpixel size in sample plane
        dx = dx_det/M
        dy = dy_det/M
    
        # coordinate system in sample plane
        # [X,Y] = meshgrid(dx*((1:1:Nx)-floor(Nx/2)-1),dy*((1:1:Ny)-floor(Ny/2)-1));
    
        # magnified coordinate system in detector plane
        # [X_det,Y_det] = meshgrid(dx_det*((1:1:Nx)-floor(Nx/2)-1),dy_det*((1:1:Ny)-floor(Ny/2)-1));
    
        # grid conversion in q-space
        dqx = 2*np.pi/(Nx*dx)
        dqy = 2*np.pi/(Ny*dy)
    
        [Qx,Qy] = np.meshgrid(dqx*(range(1,Nx+1)-np.floor(Nx/2)-1),dqy*(range(1,Ny+1)-np.floor(Ny/2)-1))
    
        ##################################
        # Prepare data from ProxToolbox:
        ##################################
        # Fresnel propagator:
        config['use_farfield_formula']=False
        config['Nx']=Nx
        config['Ny']=Ny
        config['fresnel_nr'] = 2*np.pi * (dx*dy*Nx*Ny)/(z_eff*lambd);# 1*2*pi*Nx; # Fresnel number
        config['magn']=1 # magnification factor (should this be M above, or 
        # is it okay since the demagnified pixel size is used?)

        kappa = np.sqrt(k**2-(Qx**2+Qy**2))
        config['FT_conv_kernel'] = fftshift(exp(-1j*kappa*z_eff))
        
        # problem data
        config['data'] = data['DP']
        config['rt_data'] = sqrt(data['DP'])
        config['data_zeros'] = config['data']==0
        config['norm_rt_data'] =sqrt(sum(sum(config['data'])))
        # use the abs_illumination field to represent the 
        # support constraint.
        config['abs_illumination'] = np.ones(data['DP'].shape)
        config['support_idx'] = config['abs_illumination'].nonzero()
        config['product_space_dimension'] = 1
        config['supp_phase'] = []

        # start values
        config['u_0']=ifft2(fft2(config['rt_data']*config['magn'])/config['FT_conv_kernel'])
       
        '''
        figure(1);
        subplot(2,2,1)
        imagesc(Nx*dx*linspace(0,1,Ny)*1E6,Nx*dx*linspace(0,1,Nx)*1E6,method_input.data);
        axis equal;
        axis tight;
        xlabel('x [\mum]');
        ylabel('y [\mum]');
        colormap(gray); colorbar; title('near field observation')
        subplot(2,2,2)
        imagesc(Nx*dx*linspace(0,1,Ny)*1E6,Nx*dx*linspace(0,1,Nx)*1E6,method_input.abs_illumination);
        axis equal;
        axis tight;
        xlabel('x [\mum]');
        ylabel('y [\mum]');
        colormap(gray); colorbar; title('support constraint')
    
    
        subplot(2,2,3)
        imagesc(Nx*dx*linspace(0,1,Ny)*1E6,Nx*dx*linspace(0,1,Nx)*1E6,(abs(method_input.u_0)));
        axis image; colormap(gray); colorbar;
        xlabel('x [\mum]');
        ylabel('y [\mum]');
        title('initial guess amplitude');
    
        # # pattern
        subplot(2,2,4);
        imagesc(Nx*dx*linspace(0,1,Ny)*1E6,Nx*dx*linspace(0,1,Nx)*1E6,(angle(method_input.u_0)));
        axis image; colormap(gray); colorbar;
        xlabel('x [\mum]');
        ylabel('y [\mum]');
        title('initial guess phase');
        caxis([-0.9 -0.4]);
        '''
        
    elif(experiment == 'Near_field_cell_syn'):
        
        print('Loading data file: cell256.mat')
        obj = loadmat('../InputData/Phase/cell256.mat')
        true_object = exp(1j*obj['pattern']) # phase object

        res= true_object.shape[0]
        config['Nx']=res
        config['Ny']=res
        config['use_farfield_formula']=False
        config['magn']=1 # magnification factor 
        config['fresnel_nr'] = 1*2*np.pi*config['Nx'] # Fresnel number
        [Y,X] = np.meshgrid(np.arange(-res,res,2)/(2*res),np.arange(-res,res,2)/(2*res))
        
        # Fresnel propagator:
        #   the Fresnel propogator includes all relevant physical parameters
        #   so that the numerics expert need not be concerned with this
        config['FT_conv_kernel'] = fftshift(exp(-1j*(2*np.pi*res)**2/(2*config['fresnel_nr'])* (X**2 + Y**2)))

        # simulate the hologram:
        prop_data = ifft2(config['FT_conv_kernel'] * fft2(true_object))
        config['data'] = mypoissonrnd(config['snr']*abs(prop_data)**2)/config['snr']
        config['data_zeros'] = config['data']==0
        config['rt_data'] = sqrt(config['data'])
        config['norm_rt_data'] =sqrt(sum(sum(config['data'])))
        
        # estimate the support constraint by a blurring operation.  
        #  In practice, the physicist would simply draw an envelope of
        # the hologram via visual inspection, but we do this automatically 
        # using the true object.  This leads to systemmatic errors in the 
        # reconstruction, which a methematician does not care about.
        g=gaussian(256,np.array([.1,.1]),np.array([127, 127])) #need to shift enter index by one compared to matlab
        blurr=fftshift(np.real(ifft2(fft2(obj['pattern'])*fft2(g))))
        config['support_idx'] = blurr>10
        config['indicator_ampl'] = (blurr > 10).astype(np.float32)
        config['truth'] = true_object
        config['truth_dim'] = true_object.shape
        config['norm_truth']= np.linalg.norm(true_object,'fro')
        # initial guess:
        config['u_0']=config['indicator_ampl']

        '''
        figure(1)
        subplot(2,2,1), 
        imagesc(real(true_object));
        title('True Object')
        subplot(2,2,2), 
        imagesc(method_input.rt_data); title('Near field data')
        subplot(2,2,3), 
        imagesc(method_input.indicator_ampl); title('Support constraint')
        '''
                
    elif(experiment == 'living_worm'):
        # Measuremnt of a tadpole head with the laboratory source JULIA 
    
        # hologram
        print('Loading data data_living_worm')
        
        P = loadmat('../InputData/Phase/data_living_worm.mat')
        P = P['P'][0,0]
        
        ## single image reconstruction
    
        # total number of elements
        N = P['I_exp'].size
    
        #number of pixels
        Ny = P['I_exp'].shape[0]; Nx = P['I_exp'].shape[1]
        config['Ny']=Ny; config['Nx']=Nx

        ##########################
        # Experimental parameters
        ##########################
    
        # energy in keV
        E = P['E']
    
        # wavelength [m]
        lambd = P['lambda']
        k = 2*np.pi/lambd
    
        # effective distance of detector
        z_eff = P['z_eff']
    
        # effective magnification factor
        M = P['M']
    
        #pixel size in detector plane
        dx_det = P['dx']
        dy_det = P['dx']
    
        # magnified coordinate system in detector plane
        # [X_det,Y_det] = meshgrid(dx_det*[0:1:Nx-1],dy_det*[0:1:Ny-1]);
    
        # demagnifiedpixel size in sample plane
        dx = dx_det/M
        dy = dy_det/M
    
        # coordinate system in sample plane
        # [X,Y] = meshgrid(dx*((1:1:Nx)-floor(Nx/2)-1),dy*((1:1:Ny)-floor(Ny/2)-1));
    
        # magnified coordinate system in detector plane
        # [X_det,Y_det] = meshgrid(dx_det*((1:1:Nx)-floor(Nx/2)-1),dy_det*((1:1:Ny)-floor(Ny/2)-1));
    
        # grid conversion in q-space
        dqx = 2*np.pi/(Nx*dx)
        dqy = 2*np.pi/(Ny*dy)
    
        [Qx,Qy] = np.meshgrid(dqx*(range(1,Nx+1)-np.floor(Nx/2)-1),dqy*(range(1,Ny+1)-np.floor(Ny/2)-1))
        
        ##################################
        # Prepare data from ProxToolbox:
        ##################################
        # Fresnel propagator:
        config['use_farfield_formula']=False
        config['Nx']=Nx
        config['Ny']=Ny
        config['fresnel_nr'] = P['fresnelnumber'] # Fresnel number
        config['magn']=1 # magnification factor (should this be M above, or 
        # is it okay since the demagnified pixel size is used?)

        kappa = np.sqrt(k**2-(Qx**2+Qy**2))
        config['FT_conv_kernel'] = fftshift(exp(-1j*kappa*z_eff))
        
        # problem data
        config['data'] = P['I_exp']
        config['rt_data'] = sqrt(P['I_exp'])
        config['data_zeros'] = config['data']==0
        config['norm_rt_data'] =sqrt(sum(sum(config['data'])))
        # use the abs_illumination field to represent the 
        # support constraint.
        config['abs_illumination'] = np.ones(P['I_exp'].shape)
        config['support_idx'] = config['abs_illumination'].nonzero()
        config['product_space_dimension'] = 1
        config['supp_phase'] = []

        # start values
        config['u_0']=ifft2(fft2(config['rt_data']*config['magn'])/config['FT_conv_kernel'])

        ''' 
        figure(1);
        subplot(2,2,1)
        imagesc(Nx*dx*linspace(0,1,Ny)*1E6,Nx*dx*linspace(0,1,Nx)*1E6,method_input.data);
        axis equal;
        axis tight;
        xlabel('x [\mum]');
        ylabel('y [\mum]');
        colormap(gray); colorbar; title('near field observation')
        subplot(2,2,2)
        imagesc(Nx*dx*linspace(0,1,Ny)*1E6,Nx*dx*linspace(0,1,Nx)*1E6,P.supp);
        axis equal;
        axis tight;
        xlabel('x [\mum]');
        ylabel('y [\mum]');
        colormap(gray); colorbar; title('support constraint')
    
    
        subplot(2,2,3)
        imagesc(Nx*dx*linspace(0,1,Ny)*1E6,Nx*dx*linspace(0,1,Nx)*1E6,(abs(method_input.u_0)));
        axis image; colormap(gray); colorbar;
        xlabel('x [\mum]');
        ylabel('y [\mum]');
        title('initial guess amplitude');
    
        # # pattern
        subplot(2,2,4);
        imagesc(Nx*dx*linspace(0,1,Ny)*1E6,Nx*dx*linspace(0,1,Nx)*1E6,(angle(method_input.u_0)));
        axis image; colormap(gray); colorbar;
        xlabel('x [\mum]');
        ylabel('y [\mum]');
        title('initial guess phase');
        caxis([-0.9 -0.4]);
        '''

