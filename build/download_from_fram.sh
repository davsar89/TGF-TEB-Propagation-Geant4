rm -rf ./output_fram/*
rsync -rv --min-size=20 dsarria@fram.sigma2.no://cluster/work/users/dsarria/G4_projects/TGF-TEB-propa/build/simu_data_filtered.mat ./output_fram/ 2>&1 >/dev/null
# 2>&1 >/dev/null is to silent

#cat ./output_fram/* > fused.out

#matlab -nodisplay -nodesktop -r "convert_to_mat_file"
