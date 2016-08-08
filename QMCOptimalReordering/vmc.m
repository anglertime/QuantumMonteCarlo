global constK Box datafile detRatio_dense sideLength max_dist_sqr max_Acut dim 

dim = 3;
N = 1024;
totalSteps = 100;
[ParticleConfig, OrbitalConfig, Box, sigma, constK, sideLength, max_dist_sqr, max_Acut] = read_configuration(N);
datafile = fopen(sprintf('inexact GMRES-OutputK%d.txt',constK),'wt+');
datafile1 = fopen(sprintf('inexact Large-Error-monitor-OutputK%d.txt',constK),'wt+');

% A = compute_matrix (ParticleConfig, OrbitalConfig);
A = dlmread('matrix30K.txt');
[ cutoff, ~ ] = max_min_diagonal( A, 2*N+2, 1e-3, 1e-3);
psum = sum(diag(A));

disp(sprintf('Max element in first row = %f \n', max(A(1,:))));
disp(sprintf('Max element in first column = %f \n', max(A(:,1))));
for tempCount = 1:20
    disp(sprintf('%e ', A(1, tempCount)));
end
disp(sprintf('\n'));
for tempCount = 1:20
    disp(sprintf('%e ', A(tempCount, 1)));
end
maxN = 3; 
pairHistogram(1:100) = 0; 
totalHistogram(1:maxN, 1:maxN, 1:maxN) = 0; 
totalEnergy = 0;

%Initial coordinates
initialStep = fopen('initialParticle.txt','wt+');
save_config(initialStep, N, ParticleConfig);
fclose(initialStep);

%Energy
energy = Energy(ParticleConfig, OrbitalConfig, A, constK, 1, 1);
fprintf(datafile, 'energy: %f \n', energy);

%SkHistogram
SkHistogram = CalcStructureFactor(ParticleConfig, maxN); 
for kx = 1:maxN 
   for ky = 1:maxN 
        for kz = 1:maxN 
            fprintf(datafile, 'SkHistogram(%d,%d,%d) = %f \n', kx, ky, kz, SkHistogram(kx,ky,kz));
        end 
   end 
end        

movingParticleArr(totalSteps)=0;
pass = 0;
count_each_part(N)=0;
acceptCount = 0;

nr_v(totalSteps) = 0;
nr_in(totalSteps) = 0;
test_result(totalSteps) = 0;
error(totalSteps) = 0;
error1(totalSteps) = 0;
reord = 0;
for mcStep = 1:totalSteps
    
    
   [movingParticle, move] = random_move (N, sigma);
   randNum = rand;
   movingParticleArr(mcStep) = movingParticle;

   [u, ParticleConfig1] = move_particle (A, movingParticle, move, ParticleConfig, OrbitalConfig);

   fprintf(datafile, 'At step %d \n \tmoving part %d from %s to %s \n', mcStep,...
       movingParticle, pos_to_string(ParticleConfig(movingParticle,:)), pos_to_string(ParticleConfig1(movingParticle,:)));

   %--------------  Dense Method  --------------------
   detRatio_dense = dense_detRatio(A,u,movingParticle);
   
   %--------------  Iterative Method  -------------------
%    [detRatio_iter, iter, err, err1, flag, r_nrm, innertol, inner_nr, in_iter] = ine_stgmres(1e-3, A, u, movingParticle, 1e-4);
    [detRatio_iter, iter, err, err1, flag, r_nrm, innertol, inner_nr, in_iter] = fgmres_reorder_local(1e-3, A, u, movingParticle, 1e-4);
    
    nr_v(mcStep) = iter;
    nr_in(mcStep) = in_iter;
    error(mcStep) = err;
    error1(mcStep) = err1;
    f1 = detRatio_dense*detRatio_dense;
    f2 = detRatio_iter*detRatio_iter;
    f = abs(min(f1, 1) - min(f2, 1));
    test_result(mcStep) = f;
    if (err > 3e-3)
        fprintf(datafile1, 'At step %d, error is %8.4e and f is %8.4e  \n', mcStep, err, f);
    end
   fprintf(datafile, '\tdense: %f   iterative:%f \n', detRatio_dense, detRatio_iter);
   
   detRatio = detRatio_dense;
   if (detRatio*detRatio > randNum)
       fprintf(datafile, '-----accept   %f   %f \n\n', detRatio*detRatio, randNum);
       A(movingParticle,:) = A(movingParticle,:) + u';
       ParticleConfig = ParticleConfig1;
       acceptCount = acceptCount+1;
       %reorder
       if (mod(acceptCount, floor(0.6*N)) == 0) || (min(diag(A)) < min(0.5*cutoff, 8e-5)) || (sum(diag(A)) < psum - 1)
           temp1 = min(diag(A)) - min(0.5*cutoff, 5e-4); temp2 = sum(diag(A)) - psum + 1; 
           [ cutoff, ~] = max_min_diagonal( A, 2*N+2, 1e-3, 1e-3 );
           [ A, ~, ParticleConfig, ~ ] = max_weighted_match( A, cutoff, 3, ParticleConfig );
           fprintf(datafile1, sprintf('\n At step %d, reordered the matrix. cutoff value is %f. temp12 = %f, %f \n', mcStep, cutoff, temp1, temp2));
           psum = sum(diag(A));
%            cond1(reord) = cond(A);
           reord = reord + 1;
       end
   else
       fprintf(datafile, '-----reject    %f   %f \n\n', detRatio*detRatio, randNum);
   end

    % Find if a pass is completed
    for count = 1:N
        count_each_part(count) = length(find(movingParticleArr == count));
    end

    % Check the observables every pass exactly
    if (rem(mcStep,N)== 0 && min(count_each_part) > pass)
        pass = pass + 1;
        %CalcPairCorrelation(ParticleConfig, pairHistogram);

        SkHistogram = CalcStructureFactor(ParticleConfig, maxN); 
        for kx = 1:maxN 
           for ky = 1:maxN 
                for kz = 1:maxN 
                    fprintf(datafile, 'SkHistogram(%d,%d,%d) = %f \n', kx, ky, kz, SkHistogram(kx,ky,kz));

                    totalHistogram(kx,ky,kz) = totalHistogram(kx,ky,kz) + SkHistogram(kx,ky,kz);
                end 
           end 
        end        
        energy = Energy(ParticleConfig, OrbitalConfig, A, constK, 1, 1);
        fprintf(datafile, 'energy: %f \n', energy);

        totalEnergy = totalEnergy + energy;
    end
    
    %Save configuration to check equilibrium
     if (mod(mcStep, 5000) == 0)
         Steps = fopen(sprintf('%dparticle.txt', mcStep), 'wt+');
         save_config(Steps, N, ParticleConfig);
         fclose(Steps);   
         save(sprintf('%ddata', mcStep));
         fprintf(datafile1, sprintf('\n At step %d, accept move is %d, \n Accept ratio is %6.4f%%. \n', ...
             mcStep, acceptCount, 100*acceptCount/totalSteps));
     end
end

averageEnergy = totalEnergy / pass;
fprintf(datafile, 'averageEnergy: %f \n', averageEnergy);
averageHist = totalHistogram / pass;
fprintf(datafile, 'averageHist: %f \n', averageHist);

fprintf(datafile1, '\n %d out of %d movements are accepted. Accept ratio is %6.4f%%  \n', ...
    acceptCount, totalSteps, 100*acceptCount/totalSteps);
fprintf(datafile1, '\n Finally, after %d steps, pass = %d. \n', totalSteps, pass);

save('Finaldata');          
fclose(datafile1);
fclose(datafile);