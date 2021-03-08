function B = transformTensor(A, T, condition)
%
% Tensor transformation
% 
% A = Tensor (or 3d array of tensors) in initial coordinate system (2x2 or 
%     3x3).
% T = transformation matrix, i.e. directional cosine matrix
% B = Tensor(s) in new coordinate system, same size as A.
%
if numel(size(A)) == 2
    assert(all([ismatrix(A) ~diff(size(A))]))   % Assert that A is a square matrix
    assert(sum(size(A)-size(T)) == 0)           % Assert that A and T have same size
    B = T'*A*T;
    
    if nargin>2                                 % check that condition is met.
        switch condition
            case 'posDef'                       % B must be positive definite
                assert(all(eig(B) > 0))
        end
    end
    
else
    assert(all([ismatrix(A(:, :, 1)) ~diff(size(A(:, :, 1)))]))
    assert(sum(size(A(:, :, 1))-size(T)) == 0)
    B = cellfun(@(x) T'*x*T, num2cell(A,[1 2]), 'UniformOutput', false);
    B = cat(3, B{:});
    
    if nargin > 2
        switch condition
            case 'posDef'                       % this is slow
                for n=1:size(A, 3)
                    assert(all(eig(B(:, :, n)) > 0))
                end 
        end
    end
end

end