rm -rf ./output_fram_mart/*
rsync -rv --min-size=20 dsarria@fram.sigma2.no://cluster/work/users/dsarria/G4_projects/martino_agile_rad_dist_spec/build/output_ascii/*.out ./output_fram_mart/ 2>&1 >/dev/null
# 2>&1 >/dev/null is to silent

#cat ./output_fram/* > fused.out

#matlab -nodisplay -nodesktop -r "convert_to_mat_file"
