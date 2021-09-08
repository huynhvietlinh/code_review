import csv
import pandas as pd
import numpy as np
import pybedtools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statistics import mean

def str_split_and_sort(ids):
  return ",".join([str(x) for x in sorted([int(id) for id in ids.split(',')])])

def partition_loop(in_bed_filename, debug=False):
  tmp_flank_bed_filename = "tmp_flank.bed"
  # add partition id
  loop_df = pd.read_csv(in_bed_filename, sep='\t', names=["chrom", "start", "end"])
  loop_df["partition_id"] = loop_df.index
  # generate flanks
  right_flank_df = loop_df[loop_df.end > loop_df.start]
  right_flank_df.start = right_flank_df.end - 1
  left_flank_df = loop_df[loop_df.end > loop_df.start]
  left_flank_df.end = left_flank_df.start + 1
  flank_df = pd.concat([left_flank_df, right_flank_df]).drop(columns=["partition_id"])
  # merge all flanks and write to a tmp file
  pybedtools.BedTool.from_dataframe(flank_df).sort().merge().saveas(tmp_flank_bed_filename)
  # subtract loops by flanks 
  # then merge all non-overlapped intervals into partition
  partition_df = pybedtools.BedTool.from_dataframe(loop_df)\
                                   .subtract(b=tmp_flank_bed_filename)\
                                   .sort()\
                                   .merge(c=4, o="collapse")\
                                   .to_dataframe()
  # sort the partition label
  partition_df.columns=["chrom", "start", "end", "label"]
  partition_df.label = partition_df.label.apply(str_split_and_sort)
  return partition_df
#
# Unit test
#print(partition_loop("results/sample_loop.bed", debug=True))
#print(str_split_and_sort("1,4,10,3"))
#
def comparision (experiment_name, is_count, enhancer_filename):
  med1_ctrl_filename = "../../data/2021_03_24_Med1_Peak/Med1_Ctrl.peak.bed"
  med1_kd_filename = "../../data/2021_03_24_Med1_Peak/Med1_KD.peak.bed"
  esrrb_filename = "../../data/2021_01_28_Esrrb_KD/GSM2417188_Esrrb_mm9.bed"
  all_enhancer_filename = "results/all_enhancer.bed"
  esrrb_med1_enhancer_filename = "results/esrrb_med1_enhancer.bed"
  gene_filename = "results/full_esrrb_gene_mm9.bed"

  # generate all esrrb-med1 enhancers
  all_enhancer_df = pybedtools.BedTool(enhancer_filename).sort()\
                                                .annotate(files=[med1_ctrl_filename, esrrb_filename], counts=True)\
                                                .to_dataframe()
  all_enhancer_df.to_csv(all_enhancer_filename, index=False, header=False, sep='\t')
  esrrb_med1_enhancer_df = all_enhancer_df[(all_enhancer_df.name > 0) & (all_enhancer_df.score > 0)]\
                                          [["chrom", "start", "end"]]\
                                          .to_csv(esrrb_med1_enhancer_filename, index=False, header=False, sep='\t')
  # mm9 here
  # load insulated neighborhood loops
  inact_insul_nbh_df = pd.read_csv("../../data/2021_02_10_Esrrb_neighborhood/inactive_neighborhood_mm9.csv",\
                                   skiprows=2, sep='\t', names=["refseq_id", "chrom", "start", "end"])
  #print(inact_insul_nbh_df)
  act_insul_nbh_df = pd.read_csv("../../data/2021_02_10_Esrrb_neighborhood/active_neighborhood_mm9.csv",\
                                 skiprows=2, sep='\t', names=["refseq_id", "chrom", "start", "end"])
  #print(act_insul_nbh_df)
  #tmp = act_insul_nbh_df[["chrom", "start", "end"]]
  #print(tmp.drop_duplicates())
  insul_nbh_with_gene_df = pd.concat([inact_insul_nbh_df, act_insul_nbh_df])
  #print(insul_nbh_with_gene_df)
  insl_loop_filename = "results/all_insl_loop.bed"
  insl_loop_partition_filename = "results/insl_loop_partition.bed"
  insul_nbh_with_gene_df[["chrom", "start", "end"]].drop_duplicates()\
                                                   .to_csv(insl_loop_filename, index=False, header=False, sep = '\t')
  partition_loop(insl_loop_filename, debug=True)\
                .to_csv(insl_loop_partition_filename, index=False, header=False, sep='\t')


  ## Map E-P by the new method with insulated neighborhood loops
  partition_with_enhancer_df = pybedtools.BedTool(insl_loop_partition_filename)\
                                         .annotate(files=[all_enhancer_filename, esrrb_med1_enhancer_filename], both=True)\
                                         .to_dataframe()
  partition_with_enhancer_df.columns = ["partition_chrom", "partition_start", "partition_end", "partition_id",\
                                        "all_enh_count", "all_enh_percent", "esrrb_enh_count", "esrrb_enh_percent"]
  partition_with_enhancer_df.to_csv("results/debug_partition_with_enhancer.bed", index=False, sep='\t')

  tmp = dict()
  for index, row in partition_with_enhancer_df.iterrows():
    if is_count:
      if row["partition_id"] in tmp:
        tmp[row["partition_id"]] += [row["esrrb_enh_count"], row["all_enh_count"]]
      else:
        tmp[row["partition_id"]] = [row["esrrb_enh_count"], row["all_enh_count"]]
    else:
      if row["partition_id"] in tmp:
        tmp[row["partition_id"]] += [row["esrrb_enh_percent"], row["all_enh_percent"]]
      else:
        tmp[row["partition_id"]] = [row["esrrb_enh_percent"], row["all_enh_percent"]]
 
  partition_with_gene_df = pybedtools.BedTool(gene_filename)\
                                     .intersect(b=insl_loop_partition_filename, wa=True, wb=True)\
                                     .to_dataframe()
  partition_with_gene_df.to_csv("results/debug_partition_with_gene.bed", index=False, sep='\t')
  partition_with_gene_df.columns = ["gene_chrom", "gene_start", "gene_end", "gene_symbol", "refseq_id",\
                                    "wt_max", "wt_min", "wt_exp",\
                                    "kd_max", "kd_min", "kd_avg", "log_fold",\
                                    "partition_chrom", "partition_start", "partition_end", "partition_id"]
  partition_with_gene_df["tot_esrrb"] = 0
  partition_with_gene_df["tot_enh"] = 0
  
  for i in range(len(partition_with_gene_df)) : 
    partition_with_gene_df.loc[i, "tot_esrrb"] = tmp[partition_with_gene_df.loc[i, "partition_id"]][0]
    partition_with_gene_df.loc[i, "tot_enh"] = tmp[partition_with_gene_df.loc[i, "partition_id"]][1]

  partition_with_gene_df["esrrb_enh_percent"] = partition_with_gene_df.tot_esrrb/(partition_with_gene_df.tot_enh + 1e-6)
  partition_with_gene_df.to_csv("results/debug_partition_with_gene.bed", index=False, sep='\t')

  #partition_with_gene_df = partition_with_gene_df[(partition_with_gene_df.tot_esrrb > 0) & (partition_with_gene_df.tot_enh > 0)]
  #partition_with_both_df = pd.merge(partition_with_enhancer_df, partition_with_gene_df, how="inner", on=["partition_id"])
  #partition_with_both_df["esrrb_enh_percent"] = partition_with_both_df["tot_esrrb"]/partition_with_both_df["tot_enh"]
  #partition_with_both_df = partition_with_both_df[["gene_name", "esrrb_enh_percent", "log_fold"]]\
  #                                               .sort_values(by=["esrrb_enh_percent", "log_fold"])
  #print(mean(partition_with_both_df[(partition_with_both_df.esrrb_enh_percent > 0.99) \
  #                                & (partition_with_both_df.esrrb_enh_percent < 1.01)]["log_fold"]))
  #partition_with_both_df.to_csv("results/debug_final_new_df.bed", index=False, sep='\t')

  partion_with_gene_df = partition_with_gene_df[["gene_symbol", "wt_exp", "log_fold", "esrrb_enh_percent"]].sort_values(by=["gene_symbol"])
  print(mean(partition_with_gene_df[(partition_with_gene_df.esrrb_enh_percent > 0.9) & (partition_with_gene_df.esrrb_enh_percent < 1.1)]\
             ["log_fold"]))
  partition_with_gene_df["bin"] = 0
  for i in range(len(partition_with_gene_df)):
    partition_with_gene_df.loc[i, "bin"] = round(partition_with_gene_df.loc[i, "esrrb_enh_percent"]*10)/10
  partition_with_gene_df = partition_with_gene_df.sort_values(by=["bin"])
 
  x_label = "Esrrb enhancer percentage "
  if is_count:
    experiment_name += "_peak"
    x_label += "(number of peaks)"
  else:
    experiment_name += "_bp"
    x_label += "(number of peak bps)"
  experiment_name += "_Apr_8"
  dir = "../../../../www/people/lhuynh/internal/2021/files/exp/"
  #
  partition_with_gene_df.to_csv("results/new_E_P_" + experiment_name + ".bed", index=False, sep='\t')
  partition_with_gene_df = partition_with_gene_df[partition_with_gene_df.bin < 1.01]
  #
  plt.scatter(partition_with_gene_df.esrrb_enh_percent, partition_with_gene_df.log_fold, c='b', marker='o')
  plt.axhline(y=0, linestyle='--')
  plt.xlim(-0.1,1.05)
  #plt.ylim(0,100)
  rho,p_val = stats.spearmanr(partition_with_gene_df.esrrb_enh_percent, partition_with_gene_df.log_fold)
  plt.text(0.2, 1.9, "rho=" + str(round(rho*100)/100) + ", p=" + str(p_val), fontsize = 12)
  plt.xlabel(x_label, fontsize=12)
  plt.ylabel("log2(Esrrb_KD/WT)", fontsize=12)
  
  #plt.title("Esrrb enhancer percentage vs log fold expression change of target genes",fontsize=12)
  plt.savefig(dir + "scatter_plot_" + experiment_name + ".png")
  plt.close()
  #
  print("Plot!")
  ax = sns.boxplot(x = "bin", y = "log_fold", data = partition_with_gene_df)
  ax.axhline(y=0, linestyle='--')
  ax.set(ylim=(-1,1))
  #ax = sns.swarmplot(x = "bin", y = "log_fold", data = old_esrrb_gene_df, color = ".25")
  ax.set(xlabel= x_label, ylabel= "log2_expression(Esrrb_KD/WT)")
  #ax.text(0.25, 1.3, f"Welch's t-test, p-value = {p_value}", fontsize = 9)
  ax.figure.savefig(dir + "box_plot_" + experiment_name + ".png")
  plt.close()

#comparision("our_H3K27ac_new_method", False, "../../data/2021_04_01_H3K27ac_peak/H3K27ac.our_peak.bed")
#comparision("ENCODE_mm9_new_method", False, "../../data/2021_04_01_H3K27ac_peak/H3K27ac_mm9_E14.bed")
#comparision("ENCODE_mm10_new_method", False, "../../data/2021_04_01_H3K27ac_peak/H3K27ac_liftover_mm9_E14.bed")
#comparision("our_H3K27ac_new_method", True, "../../data/2021_04_01_H3K27ac_peak/H3K27ac.our_peak.bed")
#comparision("ENCODE_mm9_new_method", True, "../../data/2021_04_01_H3K27ac_peak/H3K27ac_mm9_E14.bed")
#comparision("ENCODE_mm10_new_method", True, "../../data/2021_04_01_H3K27ac_peak/H3K27ac_liftover_mm9_E14.bed")

def split_RE_id(REs):
  if (REs != "."):
    return list(map(int, REs.split(",")))
  else:
    return []

# enhancer_df: chrom, start, end, is_target
# gene_df: chrom, start, end, log_fold
# loop_df: chrom, start, end
# mapping_mode: closest_gene, loop, hub, enhanced_hub
# experiment_name
# plot_dir
# result_dir
# parameters
def compare(enhancer_df, gene_df, loop_df, mapping_mode, experiment_name, plot_dir, result_dir, parameters):
  experiment_name = experiment_name + "_" + mapping_mode
  # map enhancers to genes
  ########################
  enhancer_df = enhancer_df.reset_index()
  enhancer_df["id"] = enhancer_df.index
  gene_df = gene_df.reset_index()
  gene_df["id"] = gene_df.index
  loop_df = loop_df.reset_index()
  loop_df["id"] = loop_df.index
  # E-P is represented by a list of enhancer index lists
  e_p = [[] for i in range(len(gene_df))]
  #
  if mapping_mode == "closest":
    gene_pos_filename = result_dir + "/gene_pos.bed"
    pybedtools.BedTool.from_dataframe(gene_df[["chrom", "start", "end", "id"]])\
                      .sort()\
                      .saveas(gene_pos_filename)
    annotated_enhancer_df = pybedtools.BedTool.from_dataframe(enhancer_df[["chrom", "start", "end", "id"]])\
                                              .sort()\
                                              .closest(gene_pos_filename)\
                                              .to_dataframe()
    annotated_enhancer_df.columns = ["enh_chrom", "enh_start", "enh_end", "enh_id",\
                                    "prom_chrom", "prom_start", "prom_end", "prom_id"]
    annotated_enhancer_df = annotated_enhancer_df[annotated_enhancer_df.prom_chrom != "."]\
                                                 .astype({"enh_id": int, "prom_id": int})
    for index, row in annotated_enhancer_df.iterrows():
      e_p[row["prom_id"]].append(row["enh_id"])
    # end mapping_mode == closest
  elif mapping_mode == "loop":
    enh_filename = result_dir + "/enh_tmp.bed" 
    prom_filename = result_dir + "/prom_tmp.bed"
    pybedtools.BedTool.from_dataframe(enhancer_df[["chrom", "start", "end", "id"]])\
                      .sort()\
                      .saveas(enh_filename)
    pybedtools.BedTool.from_dataframe(gene_df[["chrom", "start", "end", "id"]])\
                      .sort()\
                      .saveas(prom_filename)
    # Generate anchors
    window = 5000
    left_anchor_df = pd.DataFrame(data = {"chrom": loop_df.chrom,\
                                          "start": loop_df.start - window,\
                                          "end": loop_df.start + window,\
                                          "id": loop_df.index})\
                       .astype({"start": int, "end": int})
    right_anchor_df = pd.DataFrame(data = {"chrom": loop_df.chrom,\
                                           "start": loop_df.end - window,\
                                           "end": loop_df.end + window,\
                                           "id": loop_df.index})\
                        .astype({"start": int, "end": int})
    map_enh_left_anchor_df = pybedtools.BedTool.from_dataframe(left_anchor_df)\
                                               .sort()\
                                               .map(b=enh_filename, c=4, o="collapse")\
                                               .to_dataframe()\
                                               .sort_values(by=["name"])
    map_prom_left_anchor_df = pybedtools.BedTool.from_dataframe(left_anchor_df)\
                                                .sort()\
                                                .map(b=prom_filename, c=4, o="collapse")\
                                                .to_dataframe()\
                                                .sort_values(by=["name"])
    map_enh_right_anchor_df = pybedtools.BedTool.from_dataframe(right_anchor_df)\
                                               .sort()\
                                               .map(b=enh_filename, c=4, o="collapse")\
                                               .to_dataframe()\
                                               .sort_values(by=["name"])
    map_prom_right_anchor_df = pybedtools.BedTool.from_dataframe(right_anchor_df)\
                                                .sort()\
                                                .map(b=prom_filename, c=4, o="collapse")\
                                                .to_dataframe()\
                                                .sort_values(by=["name"])
    #
    for i in range(len(loop_df)):
      enhs = split_RE_id(map_enh_left_anchor_df.loc[i, "score"]) + split_RE_id(map_enh_right_anchor_df.loc[i, "score"])
      proms = split_RE_id(map_prom_left_anchor_df.loc[i, "score"]) + split_RE_id(map_prom_right_anchor_df.loc[i, "score"])
      for p in proms:
        for e in enhs:
          if (p >= len(e_p)):
            print(["Error: Promoter index out of bound", i, p, len(e_p),\
                   map_left_anchor_df.loc[i, "score"], map_right_anchor_df.loc[i, "score"]])
          e_p[p].append(e)
    for i in range(len(e_p)):
      e_p[i] = list(dict.fromkeys(e_p[i]))
    #end mapping_mode == loop
  elif mapping_mode == "hub":
    enh_filename = result_dir + "/enh_tmp.bed" 
    prom_filename = result_dir + "/prom_tmp.bed"
    pybedtools.BedTool.from_dataframe(enhancer_df[["chrom", "start", "end", "id"]])\
                      .sort()\
                      .saveas(enh_filename)
    pybedtools.BedTool.from_dataframe(gene_df[["chrom", "start", "end", "id"]])\
                      .sort()\
                      .saveas(prom_filename)
    # Generate anchors
    map_enh_loop_df = pybedtools.BedTool.from_dataframe(loop_df[["chrom", "start", "end", "id"]])\
                                        .sort()\
                                        .map(b=enh_filename, c=4, o="collapse")\
                                        .to_dataframe()\
                                        .sort_values(by=["name"])
    map_prom_loop_df = pybedtools.BedTool.from_dataframe(loop_df[["chrom", "start", "end", "id"]])\
                                         .sort()\
                                         .map(b=prom_filename, c=4, o="collapse")\
                                         .to_dataframe()\
                                         .sort_values(by=["name"])
    #
    prom_to_loop = [-1 for i in range(len(gene_df))]
    for i in range(len(loop_df)):
      proms = split_RE_id(map_prom_loop_df.loc[i, "score"])
      for p in proms:
        if prom_to_loop[p] < 0:
          prom_to_loop[p] = i
        else:
          old_i = prom_to_loop[p]
          if (loop_df.loc[i, "end"] - loop_df.loc[i, "start"]) < (loop_df.loc[old_i, "end"] - loop_df.loc[old_i, "start"]):
            prom_to_loop[p] = i
    for p in range(len(gene_df)):
      loop = prom_to_loop[p]
      if loop >= 0:
        for enh in split_RE_id(map_enh_loop_df.loc[loop, "score"]):
          e_p[p].append(enh)
  # end mapping_mode == hub 
  elif mapping_mode == "enhanced_hub":
    print("")
  else:
    print(["Error: No defined mapping mode: ", mapping_mode])
  # Evaluate E-P map
  ##################
  gene_df["target_enh_percent"] = gene_df.log_fold
  gene_df["is_non_zero"] = gene_df.log_fold
  for i in range(len(e_p)):
    all_enh_length = 0
    target_enh_length = 0
    for j in e_p[i]:
      l = enhancer_df.loc[j, "end"] - enhancer_df.loc[j, "start"]
      all_enh_length += l
      if (enhancer_df.loc[j, "is_target"] > 0):
        target_enh_length += l
    gene_df.loc[i, "target_enh_percent"] = target_enh_length/(all_enh_length + 1e-6)
    gene_df.loc[i, "is_non_zero"] = 0 if gene_df.loc[i, "target_enh_percent"] < 0.01 else 1

  x_label = "Esrrb enhancer percentage "
  # scatter plot all genes
  fig = plt.figure()
  plt.scatter(gene_df.target_enh_percent, gene_df.log_fold, c='b', marker='o')
  plt.axhline(y=0, linestyle='--')
  plt.xlim(-0.1,1.05)
  rho,p_val = stats.spearmanr(gene_df.target_enh_percent, gene_df.log_fold)
  plt.text(0.2, 1.9, "rho=" + str(round(rho*100)/100) + ", p=" + str(p_val), fontsize = 12)
  plt.xlabel(x_label, fontsize=12)
  plt.ylabel("log2(Esrrb_KD/WT)", fontsize=12)
  plt.savefig(plot_dir + "/scatter_plot_" + experiment_name + ".png")
  plt.close(fig)
  # box plot of zero vs non-zero target enhancer percentage
  fig = plt.figure()
  ax = sns.boxplot(x = "is_non_zero", y = "log_fold", data=gene_df, width=0.3)
  ax.axhline(y=0, linestyle='--')
  ax.set_xticks(range(2))
  n_zero = len(gene_df[gene_df.is_non_zero == 0])
  n_non_zero = len(gene_df[gene_df.is_non_zero > 0])
  ax.set_xticklabels(["No Esrrb \nenhancer \n(n =" + str(n_zero) + ")", "At least one \nEsrrb enhancer \n(n=" + str(n_non_zero) + ")"])
  ax.set(ylim=(-1,1))
  ax.set(xlabel= x_label, ylabel= "log2(expression(Esrrb KD)/expression(WT))")
  ax.figure.savefig(plot_dir + "/box_plot_zero_vs_non_zero_" + experiment_name + ".png")
  plt.close(fig)


# All are mm10
# Load enhancers
enhancer_pos_df = pd.read_csv("../../data/2021_04_01_H3K27ac_peak/ENCFF512HTY_mm10_E14.txt", sep='\t', names=["chrom", "start", "end"])
enhancer_df = pybedtools.BedTool.from_dataframe(enhancer_pos_df)\
                                    .sort()\
                                    .annotate(files=["results/esrrb_enhancer.mm10.bed"], counts=True)\
                                    .to_dataframe()\
                                    .rename(columns={"name": "is_target"})
# Load genes
gene_df = pd.read_csv("results/full_esrrb_gene_mm10.bed", sep='\t',\
                      names=["chrom", "start", "end", "symbol", "refseq",\
                             "ctrl_max", "ctrl_min", "ctrl_avg", "kd_max", "kd_min", "kd_avg", "log_fold"])
gene_df = gene_df.drop_duplicates(["chrom", "start", "end"], keep="last").reset_index()
# Load HiCCUP loops
hiccup_micro_c_df = pd.read_csv("../../data/2021_05_05_Micro_C_loop/hiccups_results_all/merged_loops.bedpe", skiprows=2, sep='\t',\
                                names=["chr1", "x1", "x2", "chr2", "y1", "y2", "name", "score", "strand1", "strand2", "color",\
                                       "observed", "expectedBL", "expectedDonut", "expectedH", "expectedV", "fdrBL", "fdrDonut",\
                                       "fdrH", "fdrV", "numCollapsed", "centroid1", "centroid2", "radius"])
loop_df = pd.DataFrame(data = {"chrom": "chr" + hiccup_micro_c_df.chr1,\
                               "start": round((hiccup_micro_c_df.x1 + hiccup_micro_c_df.x2)/2),\
                               "end": round((hiccup_micro_c_df.y1 + hiccup_micro_c_df.y2)/2)})\
            .astype({"start": int, "end": int})
#
cohesin_loop_df = pd.read_csv("results/cohesin_loop.mm10.bed", sep='\t', names=["chrom", "start", "end"])
#
plot_dir = "/mnt/work1/users/hoffmangroup/www/people/lhuynh/internal/2021/files/E_P_benchmark"
#mapping_modes = ["closest"]
mapping_modes = ["closest", "loop", "hub"]
for mapping_mode in mapping_modes:
  compare(enhancer_df, gene_df, loop_df, mapping_mode, "May_9", plot_dir, "results", [])
  compare(enhancer_df, gene_df, cohesin_loop_df, mapping_mode, "May_10_ChIA_PET", plot_dir, "results", [])
