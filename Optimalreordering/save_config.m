function save_config (fileHandle, N, ParticleConfig)

for count = 1:N
    for innerCount = 1:1
        fprintf(fileHandle, '%f ', ParticleConfig(count,innerCount)); 
    end
    fprintf(fileHandle, '\n'); 
end