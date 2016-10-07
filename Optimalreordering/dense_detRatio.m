function detRatio = dense_detRatio(A, u, movingParticle)

au = A'\u;
detRatio = 1 + au(movingParticle);

