// 1-ifta-01.sce
//
// Copyright 2024 Andreas Gödecke
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//        http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

load("src/utils/lib");

// Function to calculate a diffractive optical element (DOE) from a target signal image.
// 
// Parameters:
//   imagePath: The path to the target signal image.
//   iterations: The number of iterations to perform.
//   z: The propagation distance [mm].
//   dx: The pixel size in x-direction [mm].
//   dy: The pixel size in y-direction [mm].
//   lambda: The wavelength of the light source [mm].
// 
// Returns:
//   void
function ifta(imagePath, iterations, z, dx, dy, lambda)
    // Read signal
    image = imread(imagePath);          // Matrix image from values of the image
                                        // All images are converted to gray images or RGB images.
                                        // For gray images, this is a MxN unsigned char matrix;
                                        // For RGB images, this is a MxNx3 unsigned char matrix 

    image = image(:,:,1);               // imread can return a hypermatrix, so we need to select a channel
                                        // If imread reads a truecolor image, it returns
                                        // a MxNx3 hypermatrix, so for example
                                        // image(:,:,1) stands for the red channel.
                                        // For gray images, im is a MxNx1 unsigned char matrix.
    
    G0 = double(image)';                // Output of imread is in integer and for processing, double is needed
    [ny, nx] = size(G0);
    G0 = G0(:,nx:-1:1);

    // Spatial coordinates
    x = ((1:nx)-nx/2-1)*dx;
    y = ((1:ny)-ny/2-1)*dy;

    dfx = 1/(nx*dx);
    dfy = 1/(ny*dy);

    // Spatial frequencies
    fx = ((1:nx)-nx/2-1)*dfx;
    fy = ((1:ny)-ny/2-1)*dfx;

    // Fresnel phase factors
    [X, Y] = meshgrid(x, y);
    [FX, FY] = meshgrid(fx, fy);

    V = exp(%i*2*%pi/(lambda*z)*(X.*X+Y.*Y));
    W = exp(%i*2*%pi*(lambda*z)*(FX.*FX+FY.*FY));

    // Initialize complex field in the signal plane
    g0 = sqrt(G0);

    gabs = g0;
    //garg = 2*%pi*rand(g0);
    garg = zeros(g0);
    g = gabs.*exp(%i*garg);

    // Plots
    
    fig = scf(1);
    fig.color_map = gray(256);
    clf;
    
    subplot(2,2,1);
    title("Intensity in signal plane");
    subplot(2,2,2);
    title("Phase in signal plane");
    subplot(2,2,3);
    title("Intensity in DOE plane");
    subplot(2,2,4);
    title("Phase in DOE plane");

    for i = 1:iterations
        tic();
        // Fourier transformation to the DOE plane
        f = hologram_from_signal(g, V, W);
        phase = f;

        subplot(2,2,3);
        plot_intensity(f, x, y);
        subplot(2,2,4);
        plot_phase(f, x, y);
        
        // Constraints in the DOE plane
        f = set_amplitude(f, ones(f));
        
        // Fourier transformation to the signal plane
        g = signal_from_hologram(f, V, W);

        // Plot signal
        subplot(2,2,1);
        plot_intensity(g, lambda*z*fx, lambda*z*fy);
        subplot(2,2,2);
        plot_phase(g, lambda*z*fx, lambda*z*fy);

        // Constraints in the signal plane
        g = set_amplitude(g, g0);
        
        mprintf('Iteration %u: %.3fs\n', i, toc());
    end

    phasearg = atan(imag(phase), real(phase));
        
    imwrite(imnorm(255/(2*%pi)*(%pi+phasearg)), 'results/mask-n2fq-pi-256.bmp', 100);
endfunction
