function subG= subgrad_g(G_cell, k_top)
subG= cell(size(G_cell));
for i=1:length(G_cell)
    x= G_cell{i}(:);
    [~, idx]= sort(x,'descend');
    z= zeros(size(x));
    topK= min(k_top, length(x));
    z(idx(1:topK))= 1;
    subG{i}= reshape(z, size(G_cell{i}));
end
end


