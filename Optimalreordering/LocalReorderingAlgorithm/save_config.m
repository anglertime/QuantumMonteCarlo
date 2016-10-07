function save_config (fileHandle, N, ParticleConfig)

global dim

for count = 1:N
    for innerCount = 1:dim
        fprintf(fileHandle, '%4.15f ', ParticleConfig(count,innerCount)); 
    end
    fprintf(fileHandle, '\n'); 
end