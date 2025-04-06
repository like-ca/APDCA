function S_out= update_S_cell(S_in, gradS, L)
S_out= S_in;
for rr=1:length(S_in)
    S_out{rr}= S_in{rr} - (1/L)* gradS{rr};
    S_out{rr}(isnan(S_out{rr}))=0;
    S_out{rr}(abs(S_out{rr})>1e7)=0;
end
end