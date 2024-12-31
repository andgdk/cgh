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

function plot_intensity(fig, g, x, y) // BUG: Signal is always rescaled, not ideal!
    scf(fig);
    set(gca(), "auto_clear", "on");
    gabs = abs(g);
    G = gabs.*gabs;
    // grayplot(y, x, 255/max(G)*G);
    Glog = log(1E-100+G);
    Glogmin = min(Glog);
    Glogmax = max(Glog);
    grayplot(y, x, 255/(Glogmax-Glogmin)*(Glog-Glogmin));
    fig.color_map = gray(256);
    ca = gca();
    ca.isoview = "on";
endfunction