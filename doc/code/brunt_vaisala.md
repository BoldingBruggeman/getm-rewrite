title: Brunt-Vaisala frequency
Author: Karsten Bolding, Hans Burchard and Peter Holterman

Most stable results have been obtained with a weighted average
for the \(N^2\) calculation:

\begin{equation}\label{Naveraged}
\begin{array}{l}
\displaystyle
(N^2)_{i,j,k}\approx
\frac16 \Bigg(
2\frac{b_{i,j,k+1}-b_{i,j,k}}{\frac12(h^t_{i,j,k+1}+h^t_{i,j,k})}
\\ \\ \qquad\qquad\qquad
\displaystyle
+
\frac{b_{i+1,j,k+1}-b_{i+1,j,k}}{\frac12(h^t_{i+1,j,k+1}+h^t_{i+1,j,k})}
+
\frac{b_{i-1,j,k+1}-b_{i-1,j,k}}{\frac12(h^t_{i-1,j,k+1}+h^t_{i-1,j,k})}
\\ \\ \qquad\qquad\qquad
\displaystyle
+
\frac{b_{i,j+1,k+1}-b_{i,j+1,k}}{\frac12(h^t_{i,j+1,k+1}+h^t_{i,j+1,k})}
+
\frac{b_{i,j-1,k+1}-b_{i,j-1,k}}{\frac12(h^t_{i,j-1,k+1}+h^t_{i,j-1,k})}
\Bigg).
\end{array}
\end{equation}

These stability issues need to be further investigated in the future.


