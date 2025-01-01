// index.sce
//
// Copyright 2025 Andreas Gödecke
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


workingDir = getenv('cgh-01-dir', pwd());  // Optionally, set the working directory
                                           // prior to execution via console
                                           // setenv('cgh-01-dir', 'SOMEPATH')
cd(workingDir);

// Clear the console and all variables
clc;
clear;

mprintf("Generate and reconstruct CGHs with Scilab!\n\n");

// Add the source directory to the path
getd("src");

// Example: Generate a CGH and reconstruct the signal
signalImageName = "signal-256-oss.bmp";
signalImg = './data/' + signalImageName;
if isfile(signalImg) then
    try
        ifta(signalImg, 20, 1000, 5E-3, 5E-3, 632.8E-6);
        mprintf("Saved phase mask in `results/`.\n");
    catch
        mprintf("Unknown error while generating the CGH.\n");
    end
    try
        reconstruct_hologram(1000, 5E-3, 632.8E-6);
        // End of script
        mprintf("Reconstruction complete.\n");
    catch
        mprintf("Unknown error while reconstructing the signal.\n");
    end
else
    mprintf("Could not find signal target image `%s` in `data/`.\nTerminating.\n", signalImageName);
end
