// 2-doe-reconstruction.sce
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

// Function to reconstruct computer generated hologram (CGH) from a phase mask (DOE).
// 
// Parameters:
//   f: The focal length of the Fourier lens [mm].
//   p: The pixel width or height in the reconstruction plane [mm].
//   lambda: The wavelength of the light source [mm].
// 
// Returns:
//   void
function reconstruct_hologram(f, p, lambda)
    //m_block_wanted = 6;
    //n_block_wanted = 4;
    m_block_wanted = 4;     // Specify even number
    n_block_wanted = 2;     // Specify even number

    Noise = 1;              // 1 true, otherwise false
    PhaseModulation = 1;    // 1 true, otherwise false
    RowShiftFactor = 1;
    //RowShiftFactor = 0.1; // No shift applied if 1
                            // Otherwise RowShiftFactor*rand()
    envelope = 1;           // 1 true, otherwise false

    // ---------------------------------------------------------------------------
    // Read the hologram transmission function
    // ---------------------------------------------------------------------------

    // Requires atomsInstall("IPCV");
    // Read image with IPCV toolbox
    bitmapFile = "results/mask-n2fq-pi-256.bmp";
    T = double(imread(bitmapFile));

    // Scale transmission values from T_min_wanted to 1
    //T_min_wanted = 0.5;
    T_min_wanted = 10^-0.2; // Transmission = 10^(-Absorption)

    T_max_initial = max(T);
    T_min_initial = min(T);

    T = T - T_min_initial;
    T = T / T_max_initial;

    T_max = max(T);
    T = T_min_wanted + (T_max - T_min_wanted) * T;

    //T = max(T) - T;       // invert transmitting and opaque areas of the hologram
    
    // Matrix T contains the transmission values of the individual pixels in the hologram.
    // Changes could be made to the matrix at this point. For example, adding noise to simulate
    // uneven pixel transmissions. Or element-wise multiplications with a phase matrix to simulate
    // a phase-shifting effect of the pixels.

    // ---------------------------------------------------------------------------
    // Overlay the hologram transmission function with noise
    // ---------------------------------------------------------------------------

    if Noise == 1 then
        for i = 1:size(T, 'c') do
            for j = 1:size(T, 'r') do
                k = T(i,j);
                if T(i,j) < 1.2 * T_min_wanted then // Strong noise on
                                                    // absorbing and
                                                    // pixels deviating by 1.2 times the value
                    T(i,j) = T(i,j) + 0.09*rand(T(i,j), 'normal');
                    while (T(i,j) < 0 | T(i,j) > 1) do
                        T(i,j) = k + 0.09*rand(T(i,j), 'normal');
                    end
                else                                 // weak noise on
                                                    // the rest
                    T(i,j) = T(i,j) - 0.01*rand(T(i,j), 'normal');
                    while (T(i,j) > 1 | T(i,j) < 0) do
                        T(i,j) = k - 0.01*rand(T(i,j), 'normal');
                    end
                end
            end
        end
    end

    // ---------------------------------------------------------------------------
    // Multiply the hologram transmission function with a phase matrix
    // ---------------------------------------------------------------------------

    if PhaseModulation == 1 then
        // Apply to rand(x, 'normal') -> Factor:
        T = T + 0.047*rand(T, 'normal') * %i;
        // T = T + rand(T) * %i;
    end

    // ---------------------------------------------------------------------------
    // Zero-padding to increase resolution in the reconstruction plane
    // ---------------------------------------------------------------------------

    Z = zeros(T);

    T_padded = [T Z; Z Z];
    [m_padded, n_padded] = size(T_padded);
    
    // Row and column dimensions of T
    // [m, n] = size(T);

    // Alternatives:
    //
    //   T_padded = [T Z Z; Z Z Z; Z Z Z]; 
    //   [m_padded, n_padded] = size(T_padded);
    // or
    //
    //   m_padded = 512; // or another number > m
    //   n_padded = 512; // or another number > n
    //   T_padded = zeros(m_padded, n_padded);
    //   T_padded(1:m, 1:n) = T;
    ////
    // or ...

    // ---------------------------------------------------------------------------
    // Calculation of the periodic part of the reconstruction
    // ---------------------------------------------------------------------------

    //if RowShiftFactor ~=1 then
    //    m_vector = 1:1:m_padded;
    //    n_vector = 1:1:n_padded;
    //    for i=1:m do
    //        if modulo(m_vector(i), 2) == 1 then
    //            m_vector(i) = m_vector(i) + RowShiftFactor * rand();
    //        end
    //    end
    //    
    //    // fourier transform along the rows
    //    for i=1:n_padded do
    //        T_padded1(:,i) = exp(-2*%i*%pi*(0:m_padded-1)'.*.(0:m_padded-1)/m)*T_padded(:,i);
    //    end
    //    
    //    // fourier transform along the columns
    //    for j=1:m_padded do
    //        T_padded2temp = exp(-2*%i*%pi*(0:n_padded-1)'.*.(0:n_padded-1)/n_padded)*(T_padded1(j,:)).';
    //        T_padded2(j,:) = T_padded2temp.';
    //    end
    //    
    //    U_block = T_padded2;
    //else
    //    U_block = fft(T_padded);
    //end
    U_block = fft(T_padded);

    // ---------------------------------------------------------------------------
    // Periodic repetition
    // ---------------------------------------------------------------------------

    row = U_block;
    for j = 2:n_block_wanted do
        row = [row U_block];
    end

    U_p = row;
    for i = 2:m_block_wanted do
        U_p = [U_p; row];
    end

    [m_final, n_final] = size(U_p);

    m_block = m_final/m_padded;
    n_block = n_final/n_padded;

    // Matrix U_p - the complex field amplitude in the reconstruction plane without
    // the envelope - is a block matrix consisting of m_block x n_block
    // blocks U_block. The program assumes that m_block and n_block
    // are even numbers (in practice 2, 4, or 6). The following program
    // lines ensure this by extending U_p with additional blocks if necessary.

    if modulo(m_block, 2) == 1 then
        U_p = [U_p; U_p(1:m_padded,:)];
        m_final = m_final+m_padded;
    end

    if modulo(n_block, 2) == 1 then
        U_p = [U_p U_p(:,1:n_padded)];
        n_final = n_final+n_padded;
    end

    // ---------------------------------------------------------------------------
    // Calculation of the envelope
    // ---------------------------------------------------------------------------
        
    i = (0:m_final-1)-m_final/2;
    j = (0:n_final-1)-n_final/2;

    if RowShiftFactor ~=1 then
        for k=1:size(i,2) do
            if modulo(i(k), 2) == 0 then
            else
                i(k) = i(k) + RowShiftFactor * rand();
            end
        end
    end

    U_e = p.^2 * sinc(%pi*i/m_padded)'*sinc(%pi*j/n_padded);

    // ---------------------------------------------------------------------------
    // If needed: corresponding spatial frequencies
    // ---------------------------------------------------------------------------

    xi  = i/(m_padded*p);
    eta = j/(n_padded*p);

    // ---------------------------------------------------------------------------
    // If needed: corresponding spatial coordinates
    // ---------------------------------------------------------------------------

    x = lambda*f*xi;
    y = lambda*f*eta;

    // ---------------------------------------------------------------------------
    // Calculation of the field amplitude in the reconstruction plane
    // ---------------------------------------------------------------------------

    if envelope == 1 then
        U = U_e.*U_p;
    else
        U = U_p;
    end

    // ---------------------------------------------------------------------------
    // Graphical representation
    // ---------------------------------------------------------------------------

    A = abs(U);

    A(A < 5E-4) = 5E-4;

    // Calculate only a section of the amplitude
    //A = A([(size(U, 'r')/4):(size(U, 'r')*3/4-1)],[(size(U, 'c')/4):(size(U, 'c')*3/4-1)]);

    // close(winsid());
    // clf();

    //subplot(1,2,1)
    //grayplot(1:size(A, 'r'), 1:size(A, 'c'), A);
    //surf(log(A),'facecol','interp')

    //surf(log(A),log(A)+5); // color array specified
    //e=gce();
    //e.cdata_mapping='direct' // default is 'scaled' relative to the colormap
    //e.color_flag=3; // interpolated shading mode. The default is 4 ('flat' mode) for surf
    fig = scf(2);
    Sgrayplot(x, y, log(A));
    //Sgrayplot(xi, eta, log(A));

    //fig.color_map = hotcolormap(64);
    fig.color_map = gray(256);
    //addcolor(name2rgb('grey')/255);
    //colorbar(min(A)-min(A),max(A)+min(A));
    //colorbar(0,1);

    fig.children.isoview = "on";
    fig.children.x_label.text = "Spatial coordinate y'' [mm]";
    fig.children.y_label.text = "Spatial coordinate x'' [mm]";
    fig.children.sub_ticks = [0,0];
    fig.children.font_size = 4;
    fig.children.x_label.font_size = 5;
    fig.children.y_label.font_size = 5;
    fig.children.tight_limits = "on";
    fig.children.auto_scale = "off";
    //fig.children.auto_ticks = ["off","off"];
    fig.figure_position = [1000,500];
    fig.figure_size = [2000,1000];
    xlfont('SansSerif');

    Sgrayplot(x, y, log(A)); // BUG: Graph generation needs some love, without this line, the graph is empty on first execution

    //w=getdate()
    ////time = w(1)+w(2)+w(6)+w(7)+w(8)+w(9);
    //xs2png(gcf(),'Reconstruction.png')

    //subplot(1,2,2)
    //surf(A) // interpolated shading mode (color_flag == 3)

    //surf(Z,Z+0.004) // color array specified
    //e=gce();
    //e.cdata_mapping='direct' // default is 'scaled' relative to the colormap
    //e.color_flag=3; // interpolated shading mode. 
    //The default is 4 ('flat' mode) for surf

    //scf();
    //hist3d(A);

    //mesh(1:size(U, 'r'), 1:size(U, 'c'), A);

    // ---------------------------------------------------------------------------
    // Output to file
    // ---------------------------------------------------------------------------

    //fprintfMat("Reconstruction", log(A));
endfunction