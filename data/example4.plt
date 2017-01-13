set term postscript eps enhanced "Helvetica" 20
set output "example4.eps"
set size 0.875
set logscale y
set key spacing 3
set key top right
set ytics 1e+3
set yrange [5e-17:5000]
set xlabel "{/Helvetica-Italic=25 n}"
set ylabel "{/Helvetica=25 Maximum Error}"
plot "SE_Sinc_4.dat" using 1:2 with lp title "Observed error (SE)" ls 1 lw 2 pt 2, "SE_Sinc_4.dat" using 1:3 with l title "Error estimate (SE)" ls 3 lw 3, "DE_Sinc_4.dat" using 1:2 with lp title "Observed error (DE)" ls 1 lw 2 pt 4
