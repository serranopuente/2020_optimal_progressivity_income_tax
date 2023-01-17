function sub=ind2sub_homemade(sizeA, index)

% It does much the same as ind2sub, but by default returns the vector of
% subscripts that one would expect in this case.

sub=zeros(1,length(sizeA));
sub(1)=rem(index-1,sizeA(1))+1;
for j=2:length(sizeA)-1
    sub(j)=rem(ceil(index/prod(sizeA(1:j-1)))-1,sizeA(j))+1;
end
sub(length(sizeA))=ceil(index/prod(sizeA(1:length(sizeA)-1)));

end