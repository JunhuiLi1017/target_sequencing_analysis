cut -f 1 2-Targets-XGEN.C8ABDCD35B0B466EBA1ABA0950FDFB1C.bed | sort | uniq > 2_target_xgen.chr.levels.bed
vi 2_target_xgen.remove.chr.levels.bed
awk 'NR==FNR{c[$1]=$1}NR!=FNR{if(!c[$1]){print $0}}' 2_target_xgen.remove.chr.levels.bed 2-Targets-XGEN.C8ABDCD35B0B466EBA1ABA0950FDFB1C.bed | cut -f 1 | sort | uniq
awk 'NR==FNR{c[$1]=$1}NR!=FNR{if(!c[$1]){print $0}}' 2_target_xgen.remove.chr.levels.bed 2-Targets-XGEN.C8ABDCD35B0B466EBA1ABA0950FDFB1C.bed > 2_target_xgen.nochr.bed
awk 'BEGIN{FS=OFS=","}{$1="chr"$1}1' 2_target_xgen.nochr.bed > 2_target_xgen.bed
cat 2_target_xgen.bed | awk '{print $1,$2-50,$2,$1,$3,$3+50}' | xargs -n 3  | sed 's/ /\t/g' | sort -k1,1 -k2n > 2_flank_xgen.sort.bed
sort -k1,1 -k2n 2_target_xgen.bed > 2_target_xgen.sort.bed
