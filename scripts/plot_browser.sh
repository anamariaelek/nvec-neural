
# install env
mamba create -n pygenometracks -c bioconda -c conda-forge python=3.8
mamba activate pygenometracks
mamba install -c conda-forge -c bioconda pygenometracks=3.5
mamba install -c anaconda ipykernel
python -m ipykernel install --user --name=pygenometracks_env

# create config file
gtf_fn="/home/anamaria/cluster/aelek/proj/scATAC_nvec_v2/genome/Nvec_vc1.1_long.annot.gtf"
bws_atac=$( ls /home/anamaria/cluster/aelek/proj/nvec_neural/ATACSEQ/nucleosome_free_regions/bigwig/*pos.ncfree.bigwig )
bws_rna=$( ls /home/anamaria/cluster/aelek/proj/RNAseq_nvec/*bw )
fig_dir=" /home/anamaria/cluster/aelek/proj/nvec_neural/figures/raw/"
tracks=${fig_dir}"/tracks.ini"
make_tracks_file --trackFiles ${gtf_fn} ${bws_rna} ${bws_atac} -o ${tracks}

declare -A marker_genes
marker_genes["Elav"]="Nvec_vc1.1_XM_048724864.1"
marker_genes["GATA"]="Nvec_vc1.1_XM_032367288.2"
marker_genes["Islet"]="Nvec_vc1.1_XM_032362713.2"
marker_genes["LWamide"]="Nvec_vc1.1_XM_001634346.3"
marker_genes["Ncol"]="Nvec_vc1.1_XM_032363026.2"
marker_genes["FoxQ2d"]="Nvec_vc1.1_XM_032385706.2"

declare -A marker_regions
marker_regions["Elav"]="NC_064034.1:3163613-3170759"
marker_regions["GATA"]="NC_064043.1:5890694-5915400"
marker_regions["Islet"]="NC_064038.1:9286568-9339994"
marker_regions["Ncol"]="NC_064035.1:1646597-1648850"
marker_regions["FoxQ2d"]="NC_064034.1:13083406-13091596"

# Elav
"Nvec_vc1.1_XM_048724864.1"
gene="Elav"
cord="NC_064034.1:3163613-3171509"
tracks=${fig_dir}"/tracks_"${gene}".ini"
plot=${fig_dir}/pyGenomeTracks_${gene}.svg
pyGenomeTracks --tracks ${tracks} --region ${cord} --title ${gene} -o ${plot}

# Elav
"Nvec_vc1.1_XM_032367288.2"
gene="GATA"
cord="NC_064043.1:5890160-5908208"
tracks=${fig_dir}"/tracks_"${gene}".ini"
plot=${fig_dir}/pyGenomeTracks_${gene}.svg
pyGenomeTracks --tracks ${tracks} --region ${cord} --title ${gene} -o ${plot}

# Elav
"Nvec_vc1.1_XM_048724673.1"
gene="Rpamide"
cord="NC_064036.1:19018688-19025631"
tracks=${fig_dir}"/tracks_"${gene}".ini"
plot=${fig_dir}/pyGenomeTracks_${gene}.svg
pyGenomeTracks --tracks ${tracks} --region ${cord} --title ${gene} -o ${plot}

# Elav + Fox
"Nvec_vc1.1_XM_032377897.2"
gene="Pea3_Ets"
cord="NC_064035.1:18387477-18409834"
tracks=${fig_dir}"/tracks_"${gene}".ini"
plot=${fig_dir}/pyGenomeTracks_${gene}.svg
pyGenomeTracks --tracks ${tracks} --region ${cord} --title ${gene} -o ${plot}

# Fox
"Nvec_vc1.1_XM_032385706.2"
gene="FoxQ2d"
cord="NC_064034.1:13083288-13091889"
tracks=${fig_dir}"/tracks_"${gene}".ini"
plot=${fig_dir}/pyGenomeTracks_${gene}.svg
pyGenomeTracks --tracks ${tracks} --region ${cord} --title ${gene} -o ${plot}

# Fox
"Nvec_vc1.1_XM_032365096.2"
gene="Rfx"
cord="NC_064036.1:11193754-11215637"
tracks=${fig_dir}"/tracks_"${gene}".ini"
plot=${fig_dir}/pyGenomeTracks_${gene}.svg
pyGenomeTracks --tracks ${tracks} --region ${cord} --title ${gene} -o ${plot}

# Fox + Ncol

# Ncol
"Nvec_vc1.1_XM_032363026.2"
gene="Ncol3"
cord="NC_064035.1:1646594-1649400"
tracks=${fig_dir}"/tracks_"${gene}".ini"
plot=${fig_dir}/pyGenomeTracks_${gene}.svg
pyGenomeTracks --tracks ${tracks} --region ${cord} --title ${gene} -o ${plot}

# Ncol
gene="Jun"
cord="NC_064035.1:20410393-20415063"
tracks=${fig_dir}"/tracks_"${gene}".ini"
plot=${fig_dir}/pyGenomeTracks_${gene}.svg
pyGenomeTracks --tracks ${tracks} --region ${cord} --title ${gene} -o ${plot}

# Elav + Ncol
"Nvec_vc1.1_XM_032370984.2"
gene="AshA"
cord="NC_064048.1:2244989-2249425"
tracks=${fig_dir}"/tracks_"${gene}".ini"
plot=${fig_dir}/pyGenomeTracks_${gene}.svg
pyGenomeTracks --tracks ${tracks} --region ${cord} --title ${gene} -o ${plot}

# Elav + Ncol
"Nvec_vc1.1_XM_032364887.2"
gene="FoxL2"
cord="NC_064036.1:12139517-12149925"
tracks=${fig_dir}"/tracks_"${gene}".ini"
plot=${fig_dir}/pyGenomeTracks_${gene}.svg
pyGenomeTracks --tracks ${tracks} --region ${cord} --title ${gene} -o ${plot}

# all
"Nvec_vc1.1_XM_032375602.2"
gene="Insm"
cord="NC_064040.1:11172722-11175333"
tracks=${fig_dir}"/tracks_"${gene}".ini"
plot=${fig_dir}/pyGenomeTracks_${gene}.svg
pyGenomeTracks --tracks ${tracks} --region ${cord} --title ${gene} -o ${plot}
