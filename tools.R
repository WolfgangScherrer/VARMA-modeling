# ################################################################### #
# tools.R                                                             #
# (WS)                                                  juli 17, 2018 #
# ################################################################### #

#' Plot 3-dimensional arrays
#'
#' This is as simple function for plotting three dimensional structures, like
#' the impulse response function or the autocovariance function of an ARMA
#' process.
#'
#' @param x     m-by-n-by-l array, or a vector/matrix which can be coerced
#'              to an array, see below.
#' @param dim   integer vector of length 3. If this (optional) parameter is given
#'              then x is coerced to an array with \code{x = array(x, dim = dim)}.
#' @param labels.ij string which is used to label the panels.
#'              E.g.
#'
#'              \code{labels.ij = "partialdiff*y[i_*k]/partialdiff*epsilon[j_*0]"}
#'
#'              may be used for a plot of an impulse response function. The string
#'              should contain \code{"i_"}, \code{"j_"} as place holders for the
#'              indices \code{i,j} of the respective (i,j)-th sub panels.
#' @param main  a main title for the plot.
#' @param xlab  a label for the x axis.
#' @param ...   other graphical parameters. These parameters are passed on to the
#'              \code{\link[graphics]{lines}} method.
#'
#' @return x    invisibly returns the input \code{x} after being coerced to an array.
#'
#' @importFrom graphics par plot grid abline lines axis box mtext legend
#'
#' @examples
#' plot3d(rnorm(2*3*11), dim = c(2,3,11), main ='test one',
#'        labels.ij = expression(corr(x[list(i_,t+k)],x[list(j_,t)])) )
#' plot3d(array(rnorm(4*4*21), dim = c(4,4,21)), main ='test two', xlab = 'delay (s)',
#'        labels.ij = expression(epsilon[j_*0] %->% y[i_*s]) )
#' plot3d(matrix(rnorm(4*100),nrow = 4), dim = c(4,1,100), main ='test three',
#'        xlab = 'time t-1', labels.ij = expression(series~i_) )
#' @export plot3d
plot3d = function(x, dim = NULL, labels.ij = NULL,
                  main = '', xlab = 'lag (k)', ...) {
  # check input (and convert to 3D array if necessary)
  if (!is.null(dim)) {
    if (length(dim)!=3) stop('illegal (dim) parameter!')
    x = try(array(x, dim = dim))
  }
  if (inherits(x,'try-error')) stop('non compatible arguments (x,dim)!')
  if ((!is.array(x))||(length(dim(x))!=3)) stop('x must be a 3-dimensional array!')
  m = dim(x)[1]
  n = dim(x)[2]
  lag.max = dim(x)[3] -1

  if (is.null(labels.ij)) labels.ij = expression((list(i_,j_)))

  # main plot area
  opar = par(oma = c(2,1,2,0.1), mar = c(0,0,0,0),
             mgp = c(1.25, 0.15,0), tcl = -0.2,
             cex.main = 1, cex.axis = 0.6, fig = c(0, 1, 0, 1))
  plot(1:10, type = 'n', axes = FALSE)
  omd = par('omd')
  # box(which = 'plot')
  # box(which = 'outer')

  mtext(main, side=3, outer=TRUE, adj=0.5, line = 0.5, cex = 1)
  mtext(xlab, side=1, outer=TRUE, adj=0.5, line = 1.0, cex = 0.8)

  par(oma = c(0,0,0,0), mar = c(0,0,0,0)+0.25)
  for (i in (1:m)) {
    for (j in (1:n)) {
      par(fig = c(omd[1]+(j-1)*(omd[2]-omd[1])/n, omd[1]+(j)*(omd[2]-omd[1])/n,
                  omd[3] + (m-i)*(omd[4]-omd[3])/m, omd[3] + (m-i+1)*(omd[4]-omd[3])/m),
          new = TRUE)
      plot(0:lag.max, x[i,j,], ylim = range(x), type = 'n', axes = FALSE)
      grid(col = 'grey')
      abline(h = 0, col ='grey')
      lines(0:lag.max, x[i,j,],...)
      axis(2, labels = (j==1), col = 'grey40')
      axis(1, labels = (i==m), col = 'grey40')
      # the syntax for labels.ij is quite tricky
      # so we use a try() box.
      try({
        label = eval(parse(text = paste("substitute(",
                                        labels.ij,", list(i_ = i, j_ = j))",
                                        sep="")))
        legend('topright', legend = label, cex = 0.8, text.col = 'grey40',
               bty = 'n', bg = 'white', box.col = 'white', seg.len = 0)

      })
      box(col = 'grey40')
    }
  }

  par(opar)   # reset graphical parameters
  invisible(x)
}



#' Forecast Error Variance Decomposition
#'
#' Computes the Forecast Errors Variance Decomposition from a given orthogonalized
#' impulse response function.
#'
#' @param irf orthogonalized impulse response function, as computed
#'            eg. by \code{\link[MTS]{VARMAirf}} or \code{\link{SSirf}}.
#'
#'            \code{irf} may be a 3-dimensional array of dimension
#'            n-by-n-by-l or a matrix which is then coerced to such an array, see
#'            the parameter \code{dim} below.
#' @param dim integer vector of length 3. If this (optional) parameter
#'            is given then \code{irf} is coerced to an
#'            array by \code{irf = array(irf, dim=dim)}. Note that
#'            \code{dim[1]==dim[2]} must hold.
#'
#' @return
#' \item{vd}{n-by-n-by-l array which contains the forecast error variance
#' decomposition: \code{vd[i,j,h]} is the percentage of the variance of the
#' h-step ahead forecast error of the i-th component due to the j-th
#' orthogonalized shock.}
#' \item{v}{n-by-l matrix which contains the forecast error variances:
#' \code{v[i,h]} is the variance of the h-step ahead forecast error for the
#' i-th component.}
#'
#' @seealso \code{\link[MTS]{FEVdec}}.
#'
#' @examples
#' phi = matrix(c( 0.5, -0.7,  0.3, 0.75,
#'                -0.2, -0.5,  0.4,-0.90),
#' byrow = TRUE, nrow = 2)
#' theta = matrix(c(-0.1, -0.8,
#'                  0.7, -0.1),
#'                byrow = TRUE, nrow = 2)
#' k = MTS::VARMAirf(Phi = phi, Theta = theta,
#'                   lag = 24, orth = TRUE)
#' out = fevd(k$irf, dim = c(2,2,25))
#' plotfevd(out$vd)
#'
#'
#' @export
fevd = function(irf, dim = NULL) {
  # check input (and convert to 3D array if necessary)
  if (!is.null(dim)) {
    if (length(dim)!=3) stop('illegal (dim) parameter!')
    irf = try(array(irf, dim = dim))
  }
  if (inherits(irf,'try-error')) stop('non compatible arguments (irf,dim)!')
  if ((!is.array(irf))||(length(dim(irf))!=3)||(dim(irf)[1]!=dim(irf)[2]))
    stop('irf must be a (n,n,h)-dimensional array!')

  n = dim(irf)[1]
  h.max = dim(irf)[3]
  vd = apply(irf^2, MARGIN = c(1,2), FUN = cumsum)  # vd[h,i,j] = sum[l=1:h] (irf[i,j,l]^2)
  v  = apply(vd, MARGIN = c(1,2), FUN = sum)        # v[h,i] = sum[j=1:n] (vd[h,i]) = variance of the i-th component
                                                    #                                 of the h-step ahead error
  vd = vd / array(v, dim=c(h.max,n,n))
  vd = aperm(vd, c(2,3,1))
  v = t(v)
  return(list(vd = vd, v = v))
}

#' Plot a Forecast Error Variance Decomposition
#'
#' Generate a plot of a forecast error variance decomposition.
#'
#' @param vd n-by-n-by-h array, which contains the forecast error variance
#'           decomposition, as computed e.g. by \code{\link{fevd}}.
#' @param main a main title for the plot.
#' @param xlab a label for the x axis.
#' @param series.names n-dimensional vector of character strings
#'              with the names of the n times series.
#' @param col n-dimensional vector of colors.
#'
#' @return invisible copy of \code{vd}.
#'
#' @importFrom graphics rect
#' @importFrom grDevices hcl
#'
#' @examples
#' vd = array(runif(4*4*20), dim = c(4,20,4))
#' v = apply(vd, MARGIN = c(1,2), FUN = sum)
#' vd = vd / array(v, dim = c(4,20,4))
#' vd = aperm(vd, c(1,3,2))
#' plotfevd(vd)
#' # of course this sub-sample is not a valid FEVD
#' plotfevd(vd[1:2,1:2,1:10], series.names = c('x','y'), col = c('red','blue'))
#'
#' @export
plotfevd = function(vd, main = 'forecast error variance decomposition',
                    xlab = 'forecast horizon (h)', series.names, col) {
  n = dim(vd)[1]
  h.max = dim(vd)[3]

  # main plot area
  opar = par(oma = c(2, 2, 2, 0.1), mar = c(0,0,0,0),
             mgp = c(1.25, 0.15,0), tcl = -0.2, xaxs = 'i',
             cex.main = 1, cex.axis = 0.6, fig = c(0, 1, 0, 1))
  plot(1:10, type = 'n', axes = FALSE)
  omd = par('omd')
  # box(which = 'plot')
  # box(which = 'outer')

  mtext(main, side=3, outer=TRUE, adj=0.5, line = 0.5, cex = 1)
  mtext(xlab, side=1, outer=TRUE, adj=0.5, line = 1, cex = 0.8)

  par(oma = c(0,0,0,0), mar = c(0,0,0,0)+0.25)

  if (missing(col)) {
    hues = seq(15, 375, length = n + 1)
    col = hcl(h = hues, l = 65, c = 100)[1:n]
  }

  if (missing(series.names)) series.names = paste('series',1:n)

  for (i in (1:n)) {
    par(fig = c(omd[1], omd[2],
                omd[3] + (n-i)*(omd[4]-omd[3])/n, omd[3] + (n-i+1)*(omd[4]-omd[3])/n),
        new = TRUE)
    plot(c(0.5,h.max+0.5),c(0,1), type = 'n', axes = FALSE)
    grid(nx = NA, ny = NULL, col = 'gray', lty = 'solid', lwd = 1)
    for (h in (1:h.max)) {
      x = cumsum(c(0,vd[i,,h]))
      rect(h-0.4, x[-(n+1)], h+0.4, x[-1], col = col, border='black', lwd = 0.5)
    }
    axis(1, labels = (i==n))
    axis(2, labels = TRUE)
    mtext(series.names[i], side=2, outer=FALSE, adj=0.5, line = 1.1, cex = 0.8)

    box()

    legend('right', legend = series.names, fill = col, bg = 'white', cex = 0.8,
           inset = 0.01, box.col = 'white', box.lwd = 0)

    # legend('left', legend = series.names[i], inset = -0.05, horiz = TRUE, bg = col[i])
    # rect(h.max+1-0.4,0.1,h.max+1+0.4,0.9, col = col[i], border=NA)
    # text(h.max+1,0.5, series.names[i], srt = 90, cex = 0.8)
  }

  par(opar)  # reset graphical parameters
  invisible(vd)
}




#' Check the stability of a matrix polynomial
#'
#' This function checks the roots of the determinant of a polynomial
#' matrix
#' \deqn{I - a[1]z - \dots - a[p]z^p}
#' and returns \code{TRUE} if all roots are outside the unit circle. The roots
#' are determined from the eigenvalues of the companion matrix, i.e.
#' the roots are the reciprocals of the non-zero eigenvalues of the companion
#' matrix.
#'
#'
#' @param A n-by-np matrix with the polynomial coefficients \eqn{A = [a[1],\dots,a[p]]}.
#'
#' @return a scalar logical with an attribute "z" which contains the roots.
#' @export
#'
#' @seealso \code{\link[dse]{polyrootDet}}.
is.stable = function(A) {
  n = nrow(A)
  p = ncol(A)/n
  # companion matrix
  if (p>1)  {
    x = rbind(A, cbind(diag(1,n*(p-1)),matrix(0,ncol=n,nrow=(p-1)*n)))

  } else {
    x = A
  }
  z = eigen(x, only.values=TRUE)$values
  ok = (max(abs(z))<1)
  attr(ok,'z') = 1/z[abs(z)>0]
  return(ok)
}



#' Kronecker indices
#'
#' Determine the Kronecker indices, given the indices of the basis rows of the Hankel
#' matrix of the impulse response function.
#'
#' @param basis (integer) vector with the indices of the basis rows of the
#'              Hankel matrix of the impulse response.
#' @param n     (integer) dimension of the process.
#'
#' @return n-dimensional (integer) vector with the Kronecker indices.
#'
#' @seealso \code{\link[MTS]{Kronspec}}.
#'
#' @export
#'
#' @examples
#' basis2kidx(c(1,2,3), 2)
#' basis2kidx(c(1,2,4), 2)
#' \dontrun{
#' basis2kidx(c(1,2,5),2) # this is not a "nice" basis!
#' }
basis2kidx = function(basis, n) {
  # no basis elements => rank of H is zero
  if (length(basis)==0) return(integer(n))

  p = ceiling(max(basis)/n)
  m = matrix(0,nrow=n,ncol=p)
  m[basis] = 1
  kidx = apply(m,MARGIN=1,FUN=sum)
  kidx1 = apply(m,MARGIN=1,FUN=function (x) sum(cumprod(x)))
  if (any (kidx!=kidx1)) stop('this is not a nice basis!')
  return(kidx)
}


#' Construct a VARMA model from an impulse response
#'
#' The model returned is in echelon canonical form. Note that this means in
#' particular that \eqn{phi_0=theta_0}{phi[0]=theta[0]} is (in general) a
#' lower triangular matrix.
#'
#' @param k n-by-n(lag.max+1) matrix which contains the impulse response coefficients,
#'        e.g. computed by \code{\link[MTS]{VARMAirf}}.
#' @param tol tolerance used by \code{\link[base]{qr}} to estimate the rank and
#'        the basis of the Hankel matrix \eqn{H} of the impulse response function.
#'
#' @return \item{Phi}{n-by-np matrix containing the AR parameters
#' \code{Phi=[phi1,...,phip]}.}
#' \item{Theta}{n-by-np matrix containing the MA parameters
#' \code{Theta=[theta1,...,thetap]}.}
#' \item{Phi0}{n-by-n matrix, lag zero coefficient \code{Phi0=phi0=theta0}.}
#' \item{kidx}{vector of Kronecker indices of the impulse response function.}
#' \item{Hrank}{estimated rank of the Hankel matrix, as computed by
#' \code{\link[base]{qr}}.}
#' \item{Hpivot}{vector of pivot elements, returned by  \code{\link[base]{qr}}. Note that
#' the first \code{Hrank} elements of this vector are the indices of the basis rows
#' of the Hankel matrix.}
#'
#' @seealso \code{\link[MTS]{Kronspec}}, \code{\link[MTS]{Kronid}},
#'          \code{\link[MTS]{Kronfit}} and \code{\link{impresp2SS}}.
#'
#' @importFrom stats lsfit
#'
#' @export
impresp2PhiTheta = function(k, tol = 1e-8) {
  n = nrow(k)
  lag.max = ncol(k)/n
  f = ceiling(lag.max/2) # f <=> future
  p = lag.max - f        # p <=> past
  lag.max = lag.max -1

  if (any(k[,1:n,drop=FALSE] != diag(1,n)))
    stop('lag zero coefficient must be the identity matrix!')

  # construct Hankel matrix of impulse response coefficients
  H = matrix(0, nrow = f*n, ncol = p*n)
  for (i in (1:f)) H[((i-1)*n+1):(i*n),] = k[,(i*n+1):((i+p)*n)]
  # print(H)
  # N.B. consider transposed Hankel matrix
  H = t(H)
  # QR decomposition
  qr.H = qr(H, LAPACK = FALSE, tol=tol)

  # cat('H ****************\n')
  # print(str(H))
  # print(svd(H)$d)

  # rank of H is zero!
  if (qr.H$rank ==0) {
    kidx = integer(n)
    return(list(Phi = matrix(0,nrow=n,ncol=0),
                Theta = matrix(0,nrow=n,ncol=0), Phi0 = diag(1,n),
                kidx = kidx, Hrank = qr.H$rank, Hpivot = qr.H$pivot))
  }

  # index of basis rows
  basis = qr.H$pivot[1:qr.H$rank]
  # print(str(H))
  # print(basis)
  kidx = basis2kidx(basis,n)
  # cat('kidx:',kidx,'\n')
  fpar = MTS::Kronspec(kidx, output = FALSE)
  #  fpar = kidx2freepar(kidx, model = 'varma')
  #  print(fpar)
  p = max(kidx)

  if ((ncol(H) %/% n)<(p+1))
    stop('Hankel matrix is too small, try to increase number of lags!')

  # Toeplitz matrix of IRF coefficients
  Tk = matrix(0, nrow=n*p,ncol=n*p)
  for (i in (1:p)) Tk[((i-1)*n+1):(i*n),((i-1)*n+1):(p*n)] = k[,(n+1):((p-i+2)*n)]
  #  print(Tk)
  #  a = fpar[,1:(n*(p+1)),drop = FALSE]
  #  b = fpar[,(n*(p+1)+1):(n*(2*p+1)),drop=FALSE]
  a = fpar$PhiID
  a[a==2] = NA
  b = matrix(0,nrow=n, ncol=p*n)
  for (r in (1:n)) {
    # Achtung auf Anordnung der Zeilen H[i] von H, bzw. der AR Koeffizienten a[i]
    # Koeffizienten von A werden im wesentlichen bestimmt durch
    # a[p] H[1] + a[p-1] H[2] + ... + a[1] H[p] = H[p+1]
    ind = matrix(1:(n*(kidx[r]+1)),nrow=n,ncol=(kidx[r]+1))
    ind = ind[,(kidx[r]+1):1]
    # print(ind[])
    i = which(is.na(a[r,]))
    if (length(i)>0) {
      j = ind[i]
      # print(j)
      # print(ind[r])
      # print(str(H))
      # print(lsfit(H[,j,drop=FALSE],H[,ind[r]],intercept=FALSE)$coef)
      a[r,i] = lsfit(H[,j,drop=FALSE],H[,ind[r]],intercept=FALSE)$coef
    }
    a[r,r] = -1
    if (kidx[r]>0) b[r,1:(kidx[r]*n)] =
      a[r,1:(n*kidx[r])] %*% Tk[1:(n*kidx[r]),1:(kidx[r]*n)]
  }

  a0 = -a[,1:n,drop=FALSE]
  a = a[,(n+1):(n*(p+1)),drop=FALSE]
  b = (b+a)

  if (all(a0==diag(n))) a0 = NULL

  return(list(Phi = a, Theta = b, Phi0 = a0,
              kidx = kidx, Hrank = qr.H$rank, Hpivot = qr.H$pivot))
}



#' Construct a dse::ARMA model
#'
#' This function constructs an \code{\link[dse]{ARMA}} object for given
#' AR/MA parameter matrices (in the style of the \code{MTS} package).
#'
#' The \code{dse} Package uses a simple strategy to impose restrictions
#' on the parameters of a models. It simply treats all coefficients equal
#' to one or zero as fixed and the others as "free". However for a model
#' in \emph{echelon form} in addition \eqn{a_0=b_0}{a[0]=b[0]} must hold. Ignoring this constraint
#' has two effects. First \code{\link[dse]{informationTests}} does not use the
#' correct number of free parameters and hence it does not return the correct values for
#' information criteria like AIC. Second, if the model is feed into
#' \code{\link[dse]{estMaxLik}} then the estimated model will in general not be
#' in echelon canonical form.
#'
#' In order to circumvent the first problem this function offers the switch \code{fix}.
#' For \code{fix=TRUE} the function calls \code{\link[dse]{fixConstants}} and sets
#' all zero/one entries as fixed and in addition all entries of \eqn{b_0}{b[0]}.
#' By this trick \code{\link[dse]{informationTests}} then uses the right number
#' of free parameters corresponding to a model in echelon form. However
#' this trick does not solve the second problem.
#'
#'
#' @param Phi   n-by-np matrix with the AR parameters \code{Phi=[phi1,...,phip]}.
#' @param Theta n-by-nq matrix containing the MA parameters
#'              \code{Theta=[theta1,...,thetaq]}.
#' @param Phi0  n-by-n matrix containing the lag zero coefficient matrix
#'              \code{phi0 = theta0}. In the default case \code{Phi0=NULL}
#'              the n-dimensional identity matrix is used.
#' @param fix   see details.
#' @param output.names vector of strings with the respective names for the
#'                     n components of the ARMA process.
#'
#' @return \code{\link[dse]{ARMA}} object.
#'
#' @seealso \code{\link{ARMA2PhiTheta}}.
#'
#' @examples
#' phi = matrix(c(-0.5, -0.2, -0.3, -0.05, -0.2, -0.5, -0.1, -0.30),
#' byrow = TRUE, nrow = 2)
#' theta = matrix(c(-0.2,  0.0, -0.1, -0.3),
#'                byrow = TRUE, nrow = 2)
#' phi0 = matrix(c(1.0, 0, 0.4, 1),
#'               byrow = TRUE, nrow = 2)
#' arma = PhiTheta2ARMA(solve(phi0,phi), solve(phi0,theta))
#' arma = PhiTheta2ARMA(phi, theta, phi0)
#' m = ARMA2PhiTheta(arma, normalizePhi0 = FALSE)
#' all.equal(cbind(phi0,phi,theta),cbind(m$Phi0,m$Phi,m$Theta))
#'
#'
#' @export
PhiTheta2ARMA = function(Phi, Theta, Phi0 = NULL, fix = FALSE,
                         output.names = paste('y',1:nrow(Phi),sep='')) {
  n = nrow(Phi)
  p = ncol(Phi)/n
  q = ncol(Theta)/n
  if (is.null(Phi0)) {
    a0 = diag(1,n)
  } else {
    a0 = Phi0
  }
  A = cbind(a0,-Phi)
  dim(A) = c(n,n,p+1)
  A = aperm(A, c(3,1,2))
  B = cbind(a0,-Theta)
  dim(B) = c(n,n,q+1)
  B = aperm(B, c(3,1,2))
  arma = dse::ARMA(A = A, B = B, C=NULL, output.names = output.names)
  if (fix) {
    A0 = (A == 0) | (A == 1)
    B0 = (B == 0) | (B == 1)
    B0[1,,] = TRUE
    arma = dse::fixConstants(arma, constants = list(A = A0, B = B0))
  }
  return(arma)
}


#' ARMA2PhiTheta
#'
#' Extracts the AR/MA coefficients of a \code{dse::ARMA} object and returns these
#' coefficient matrices in the style of the \code{MTS} package.
#'
#' @param arma \code{dse::ARMA} object.
#' @param normalizePhi0 if \code{TRUE} then the lag zero coefficient matrix is normalized to the
#'        identity matrix and the other coefficients are correspondingly transformed.
#'
#' @return \item{Phi}{n-by-np matrix containing the AR parameters
#' \code{Phi=[phi1,...,phip]}.}
#' \item{Theta}{n-by-nq matrix containing the MA parameters
#' \code{Theta=[theta1,...,thetaq]}.}
#' \item{Phi0}{n-by-n matrix containing the lag zero coefficient matrix \code{phi0=theta0}.}
#'
#' @seealso \code{\link{PhiTheta2ARMA}}.
#'
#' @examples
#' A = array(c(1, .5, .3, 0.4, .2, .1, 0, .2, .05, 1, .5, .3) ,c(3,2,2))
#' B = array(c(1, .2, 0.4, .1, 0, 0, 1, .3), c(2,2,2))
#' arma = dse::ARMA(A = A, B = B, C=NULL,
#'                  output.names = paste("y", 1:2, sep = ""))
#' ARMA2PhiTheta(arma)
#' m = ARMA2PhiTheta(arma, normalizePhi0 = FALSE)
#' arma.test = PhiTheta2ARMA(m$Phi, m$Theta, m$Phi0)
#' all.equal(arma, arma.test)
#'
#' @export
ARMA2PhiTheta = function(arma, normalizePhi0 = TRUE) {
  if (dse::is.TSestModel(arma)) arma = dse::TSmodel(arma)
  if(!dse::is.ARMA(arma)) stop('input argument must be an dse::ARMA object')
  A = arma$A
  B = arma$B
  if (any(A[1,,] != B[1,,])) stop('the lag zero coefficients a0, b0 must be identical')
  a0 = A[1,,]
  A = A[-1,,,drop=FALSE]
  A = aperm(A,c(2,3,1))
  dim(A) = c(dim(A)[1], dim(A)[2]*dim(A)[3])
  B = B[-1,,,drop=FALSE]
  B = aperm(B,c(2,3,1))
  dim(B) = c(dim(B)[1], dim(B)[2]*dim(B)[3])
  if (normalizePhi0) {
    return(list(Phi = -solve(a0,A), Theta = -solve(a0,B), Phi0 = diag(1,dim(A)[1])))
  } else {
    return(list(Phi = -A, Theta = -B, Phi0 = a0))
  }
}

#' Lyapunov equation
#'
#' Solve the Lyapunov equation:
#' \deqn{P = A P A' + Q}
#' It is assumed that Q is symmetric! Note: the stability of \eqn{A}
#' is not checked.
#'
#' The procedure uses the Schur decomposition method as in \cite{Kitagawa,
#' An Algorithm for Solving the Matrix Equation \eqn{X = F X F' + S},
#' International Journal of Control, Volume 25, Number 5, pages 745--753
#' (1977)}.
#'
#' Column-by-column solution method as suggested in
#' \cite{Hammarling, Numerical Solution of the Stable, Non-Negative
#' Definite Lyapunov Equation, IMA Journal of Numerical Analysis, Volume
#' 2, pages 303--323 (1982)}.
#'
#' @param A matrix
#' @param Q matrix
#'
#' @return P solution of Lyapunov equation
#' @export
#'
#' @seealso \code{\link[QZ]{qz.zgees}} is used for computing the Schur decomposition
#' of \eqn{A}.
#'
#' @examples A = matrix(rnorm(4),nrow=2,ncol=2)
#' Q = matrix(rnorm(4),nrow=2,ncol=2)
#' Q = Q %*% t(Q)
#' P = lyap(A,Q)
lyap = function(A,Q) {
  n = dim(A)[1]
  if (any(n!=c(dim(A),dim(Q)))) stop('A, Q must be square matrices')

  # transform to complex upper diagonal matrix!
  schur = QZ::qz.zgees(matrix(as.complex(A),nrow=n,ncol=n))
  # print(schur$T)
  # print(QZ::H(schur$VS) %*% A %*% schur$VS - schur$T)
  A = schur$T
  Q = QZ::H(schur$VS) %*% Q %*% schur$VS

  # print(Q)
  if (n>1) {
    for (i in (n:2)) {
      Q[i,i] = Q[i,i]/(1-abs(A[i,i])^2)
      Q[1:(i-1),1:i] = Q[1:(i-1),1:i] + (Q[i,i]*A[(1:(i-1)),i,drop=FALSE]) %*% QZ::H(A[1:i,i,drop=FALSE])
      Q[1:(i-1),i] = solve(diag(i-1)-Conj(A[i,i])*A[1:(i-1),1:(i-1),drop=FALSE],
                           Q[1:(i-1),i,drop=FALSE])
      Q[i,1:(i-1)] = Conj(Q[1:(i-1),i])
      h = A[1:(i-1),1:(i-1),drop=FALSE] %*% Q[1:(i-1),i,drop=FALSE] %*% QZ::H(A[1:(i-1),i,drop=FALSE])
      Q[1:(i-1),1:(i-1)] = Q[1:(i-1),1:(i-1)] + h + QZ::H(h)
      # print(Q)
    }
  }
  Q[1,1] = Q[1,1] = Q[1,1]/(1-abs(A[1,1])^2)
  # print(Q)

  Q = Re(schur$VS %*% Q %*% QZ::H(schur$VS))
  Q = (Q+t(Q))/2
  return(Q)
}


#' Impulse response function of a state space model
#'
#' This function computes the impulse response function and the
#' orthogonalized impulse response function of a state space model
#' in innovation form. The syntax and the return value is analogous to
#' \code{\link[MTS]{VARMAirf}}.
#'
#' The orthognalized innovations are defined by the \emph{symmetric} square root of the
#' covariance matrix \code{Sigma}. Note that the stability of the state space model
#' is not checked.
#'
#' @param ss \code{\link[dse]{SS}} object, which represents an innovation
#'           form state space model.
#' @param Sigma covariance of the innovations (n-by-n symmetric, positive
#'              definite matrix). If \code{NULL} then the \eqn{n}-dimensional identity
#'              matrix is used.
#' @param lag.max maximum lag (integer)
#' @param orth if \code{TRUE} then the orthogonalized impulse response function
#'             is computed.
#'
#' @return \item{psi}{n-by n(lag.max+1) matrix with the impulse response
#' coefficients.}
#' \item{irf}{nn-by-(lag.max+1) matrix with the orthogonalized
#' impulse response coefficients. If \code{!orth} then
#' \code{irf} is just a "reshaped" version of \code{psi}.
#' Note that \code{psi} and \code{irf} have different dimensions.}
#'
#' @export
SSirf = function(ss, Sigma = NULL, lag.max = 12L, orth = TRUE) {
  if (dse::is.TSestModel(ss)) ss = dse::TSmodel(ss) # extract model
  if (!dse::is.innov.SS(ss))
    stop('this is only implemented for innovation form state space models!')

  n = ncol(ss$K)

  psi = matrix(0,nrow = n, ncol = (lag.max+1)*n)
  psi[,1:n] = diag(1,n)
  FK = ss$K

  for (i in (1:lag.max)) {
    psi[,(i*n+1):((i+1)*n)] = ss$H %*% FK
    FK =  ss$F %*% FK
  }

  irf = psi
  if (orth) {
    if (is.null(Sigma)) Sigma = diag(1,n)
    evd = eigen(Sigma)
    Sigma2 = evd$vectors %*% diag(sqrt(evd$values)) %*% t(evd$vectors)

    for (i in (0:lag.max)) irf[,(i*n+1):((i+1)*n)] = psi[,(i*n+1):((i+1)*n)] %*% Sigma2
  }
  dim(irf) = c(n^2,lag.max+1)

  return(list(psi = psi, irf = irf))
}

#' Autocovariance and autocorrelation function of a state space model
#'
#' This function computes the autocovariance and autocorrelation function
#' of a stationary process represented by a state space model in innovation form.
#' The syntax and the return value is analogous to
#' \code{\link[MTS]{VARMAcov}}.
#'
#' Note that the stability of the state space model is not checked.
#'
#' @param ss \code{\link[dse]{SS}} object, which represents an innovation
#'           form state space model.
#' @param Sigma covariance of the innovations (n-by-n symmetric,
#'           positive definite matrix). If \code{NULL} then the
#'           \eqn{n}-dimensional identity matrix is used.
#' @param lag.max maximum lag (integer)
#'
#' @return \item{autocov}{n-by-n(lag.max+1) matrix with the autocovariances.}
#' \item{ccm}{n-by-n(lag.max+1) matrix with the autocorrelations.}
#'
#' @export
SScov = function(ss, Sigma = NULL, lag.max = 12L) {
  if (dse::is.TSestModel(ss)) ss = dse::TSmodel(ss) # extract model
  if (!dse::is.innov.SS(ss))
    stop('this is only implemented for innovation form state space models!')

  n = ncol(ss$K)

  if (is.null(Sigma)) Sigma = diag(1,n)
  P = lyap(ss$F, ss$K %*% Sigma %*% t(ss$K))
  M = ss$F %*% P %*% t(ss$H) + ss$K %*% Sigma

  autocov = matrix(0, nrow = n, ncol = n*(lag.max+1))
  autocov[,1:n] = ss$H %*% P %*% t(ss$H) + Sigma
  sd = sqrt(diag(autocov[,1:n]))
  sd = outer(sd,sd,'*')
  ccm = autocov
  ccm[,1:n] = ccm[,1:n] /sd

  for (i in (1:lag.max)) {
    autocov[,(i*n+1):((i+1)*n)] = ss$H %*% M
    ccm[,(i*n+1):((i+1)*n)] = autocov[,(i*n+1):((i+1)*n)]/sd
    M =  ss$F %*% M
  }

  return(list(autocov = autocov, ccm = ccm))
}


#' Construct a state space model from an impulse response.
#'
#' This function implements the Ho-Kalman realization algorithm which determines a state
#' space realization for a given impulse response function.
#'
#' For the case \code{type=="echelon"} a state space system in echelon canonical form
#' is returned and for \code{type=="balanced"} a kind of balanced realisation. The core
#' step of this algorithm is to determine a basis for the row space of
#' the Hankel matrix \eqn{H} of the impulse response coefficients.
#'
#' In the first case this is done via a QR decomposition of the transposed Hankel
#' matrix with the R function \code{\link{qr}} using the tolerance parameter \code{tol}.
#' If the optional parameter \code{s} is missing or \code{NULL} then the rank of the
#' Hankel matrix and hence the state space dimension is determined from this
#' QR decomposition. If \code{s} is given then the algorithm constructs a state
#' space model with this specified state space dimension. However, this option should
#' be used with care. There is no guarantee that the so constructed model is
#' a reasonable approximation for the true model.
#'
#' For the "balanced" case an SVD decomposition of the Hankel matrix is used. For missing
#' \code{s} the rank of the Hankel matrix is determined from the singular values of
#' the matrix using \code{tol} as a threshold, i.e. the rank is estimated as the number of
#' singular values which are larger than or equal to \code{tol} times the
#' largest singular value. If \code{s} is specified then the first \code{s} right
#' singular vectors are used as a basis for row space of the Hankel matrix and a state space
#' model with state space dimension \code{s} is returned.
#'
#' @param k n-by-n(lag.max+1) matrix which contains the impulse response coefficients,
#'        e.g. computed by \code{\link[MTS]{VARMAirf}}.
#' @param type (string) determines the "realization type".
#' @param tol tolerance parameter, see details.
#' @param s desired state dimension, see details.
#'
#' @return \item{ss}{a \code{\link[dse]{SS}} object which represents the
#'         constructed innovation form state space model.}
#' \item{Hsv}{for \code{type=="balanced"}: a vector with the singular values
#' of the Hankel matrix of the impulse response.}
#' \item{kidx}{for \code{type=="echelon"}: vector of Kronecker indices of
#' the impulse response function.}
#' \item{Hrank}{for \code{type=="echelon"}: estimated rank of the Hankel matrix,
#' as computed by \code{\link[base]{qr}}.}
#' \item{Hpivot}{for \code{type=="echelon"}: vector of pivot elements, returned by
#' \code{\link[base]{qr}}. Note that the first \code{Hrank} elements of this vector
#' are the indices of the basis rows of the Hankel matrix.}
#'
#' @seealso \code{\link{impresp2PhiTheta}}.
#'
#' @export
impresp2SS = function(k, type = c('echelon','balanced'), tol = 1e-8, s) {
  type = match.arg(type)
  n = nrow(k)
  lag.max = ncol(k)/n
  f = ceiling(lag.max/2) # f <=> future
  p = lag.max - f        # p <=> past
  lag.max = lag.max -1

  if (any(k[,1:n,drop=FALSE] != diag(1,n))) stop('lag zero coefficient must be the identity matrix!')

  if (missing(s)) s = NULL
  if (!is.null(s)) {
    s = as.integer(s[1])
    if ((min((f-1),p)*n)<s) stop('Hankel matrix is too small for the desired state dimension!')
  }

  # construct Hankel matrix of impulse response coefficients
  H = matrix(0, nrow = f*n, ncol = p*n)
  for (i in (1:f)) H[((i-1)*n+1):(i*n),] = k[,(i*n+1):((i+p)*n)]

  if (type=='balanced') {
    H1 = H[(n+1):(f*n),,drop=FALSE]
    H = H[1:((f-1)*n),,drop=FALSE]
    #  print(H)

    # "balanced" realization via SVD
    svd.H = svd(H)
    if (is.null(s)) {
      # determine state dimension from singular values
      s = ifelse(svd.H$d[1]>.Machine$double.eps, sum(svd.H$d >= (tol*svd.H$d[1])),0)
    }
    if (s>0) {
      sv2 =  sqrt(svd.H$d[1:s])
      SH = matrix(sv2,nrow=s,ncol=n*p,byrow=FALSE) * t(svd.H$v[,1:s,drop=FALSE])
      SH1 = (matrix(1/sv2,nrow=s,ncol=n*(f-1),byrow=FALSE)*t(svd.H$u[,1:s,drop=FALSE])) %*% H1
      AC = t(lsfit(t(SH),t(rbind(SH1,H[1:n,,drop=FALSE])),intercept=FALSE)$coef)
      B = SH[,1:n,drop=FALSE]
      ii = 1:s
    } else {
      AC = matrix(0,nrow=n,ncol=0)
      B =  matrix(0,nrow=0,ncol=n)
      ii = integer(0)
    }
    return(list(ss = dse::SS(F = AC[ii,,drop=FALSE], K = B,
                             H = AC[(s+1):(s+n),,drop=FALSE]),
                Hsv = svd.H$d))
  }

  # "Echelon form" realization via QR
  qr.H = qr(t(H[1:((f-1)*n),,drop=FALSE]), LAPACK = FALSE, tol=tol)
  if (is.null(s)) s = qr.H$rank
  if (s>0) {
    # index of basis rows
    basis = qr.H$pivot[1:s]

    kidx = basis2kidx(basis,n) # Kronecker indices

    AC = t(lsfit(t(H[basis,,drop=FALSE]),t(H[c(basis+n,1:n),,drop=FALSE]),intercept=FALSE)$coef)
    dimnames(AC) = NULL

    # impose exact echelon form structure
    ind = c(basis+n,1:n)
    for (i in (1:(s+n))) {
      j = which(ind[i]==basis)
      if (length(j)==1) {
        AC[i,] = 0
        AC[i,j] = 1
      } else {
        j = which( ind[i] < basis )
        if (length(j)>0) AC[i,j] = 0
      }
    }

    B = H[basis,1:n,drop=FALSE]
    ii = 1:s
  } else {
    kidx = integer(n)
    AC = matrix(0,nrow=n,ncol=0)
    B = matrix(0,nrow=0,ncol=n)
    ii = integer(n)
  }
  return(list(ss = dse::SS(F = AC[ii,,drop=FALSE], K = B,
                           H = AC[(s+1):(s+n),,drop=FALSE]),
              kidx = kidx, Hrank = qr.H$rank, Hpivot = qr.H$pivot))
}


