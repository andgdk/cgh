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

function plot_phase(fig, g, x, y)
    scf(fig);
    set(gca(), "auto_clear", "on");
    garg = atan(imag(g), real(g));
    grayplot(y, x, 255/(2*%pi)*(%pi+garg))
    fig.color_map = gray(256);
    ca = gca();
    ca.isoview = "on";
endfunction