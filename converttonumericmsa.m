function seqs = converttonumericmsa(msa)
seqs =zeros(length(msa), length(msa{1}));
for i=1:length(msa)
  seqs(i,:) = my_aa2int(msa{i});
  %seqs(i,(seqs(i,:)==25))=21;
end;
