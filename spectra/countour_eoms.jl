using Plots
using LaTeXStrings

f(a, b) = sin(a - b) + cos(a - b)/((cos(a+b) + sin(a+b)))

a_range = 0:0.01:π/4
b_range = 0:0.01:π/4

zlim= (0, 0.25)
plot(a_range, b_range, f, st = :surface, xlabel = L"\alpha", ylabel = L"\beta", zlabel = L"f(\alpha, \beta)", interactivity = true, camera = (30, 30), zlims = zlim)
