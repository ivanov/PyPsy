
xidxs=numpy.load('xidxs.npy')
yidxs=numpy.load('yidxs.npy')[90:]
log_param_err = numpy.load('ssq_params200.npy')
contourf( xidxs, yidxs, log_param_err[:,90:].T )
#colorbar()
xlabel('p[0] (threshold)')
ylabel('p[1] (slope or JND)')
title('log(sse of optimized params)')

plot( 0.058, 0.015, '*')

