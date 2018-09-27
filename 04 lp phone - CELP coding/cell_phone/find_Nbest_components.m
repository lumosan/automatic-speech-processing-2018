function [gains, indices] = find_Nbest_components(signal, ...
    codebook_vectors, N)

% function [gains, indices] = find_Nbest_components(signal, ...
%                             codebook_vectors, codebook_norms , N)
%
% This function finds the N best components of signal from the 
% vectors in codebook_vectors, so that the residual error:
%    error=signal- codebook_vectors(indices)*gains 
% is minimized.
% Components are found one-by-one using a greedy algorithm. 
% When components in codebook_vectors are not orthogonal, 
% the search is therefore suboptimal. 

[M, L] = size(codebook_vectors);

for j = 1:L
        codebook_norms(j) = norm(codebook_vectors(:,j));
end

gains = zeros(N,1);
indices = ones(N,1);

for k = 1:N
    max_norm = 0;
    for j = 1:L
        beta = codebook_vectors(:,j)'*signal;
        if codebook_norms(j) ~= 0
            component_norm = abs(beta)/codebook_norms(j);
        else
            component_norm = 0;
        end
        if component_norm  > max_norm
            gains(k) = beta/(codebook_norms(j)^2);
            indices(k) = j;
            max_norm = component_norm ;
        end
    end
    signal = signal - gains(k)*codebook_vectors(:,indices(k));
end
