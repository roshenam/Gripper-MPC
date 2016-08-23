load('Inv_Set_Data.mat')
idx = zeros(size(DATA));
                                    
%%

counter = 0;
N30data = {};
N30params = {};

for i=1:length(DATA)
    idx(i) = DATA(i).failure;
    if isa(DATA(i).data,'struct')
        if isa(DATA(i).data.final,'struct')
            if DATA(i).params.N == 30
                counter = counter + 1;
                N30data{counter} = DATA(i).data;
                N30params{counter} = DATA(i).params;
            end
        end
    end    
end

%% 

for i=1:length(N10params)
    if N10p{i}