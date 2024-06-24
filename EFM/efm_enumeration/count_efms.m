model_path='./matlab_models/';
filenames=dir(model_path);

counts=cell(18,2);
for i=3:20
    cur_filename=filenames(i).name;
    load([model_path,cur_filename])
    reaction_ids=cellstr(reaction_ids);
    metabolite_ids=cellstr(metabolite_ids);
    stru=struct();
    stru.stoich=S;
    stru.reversibilities=reversibilities;
    stru.reactionNames=(reaction_ids)';
    stru.metaboliteNames=(metabolite_ids)';
    
    disp(['Counting modes for ',cur_filename])
    disp('.................................................')
    disp('\n\n\n\n')
    opts=CreateFluxModeOpts('count-only',true,'arithmetic', 'fractional');
    %opts=CreateFluxModeOpts('count-only',false,'arithmetic', 'fractional','enforce','enforce.txt');
    mnet = CalculateFluxModes(stru,opts);
    counts{i,1}=cur_filename;
    counts{i,2}=mnet;
end

%%
save('./efm_counts/counts.mat','counts')
%%
table=cell2table(counts);
writetable(table,"./efm_counts/efm_counts_unfiltered.csv")