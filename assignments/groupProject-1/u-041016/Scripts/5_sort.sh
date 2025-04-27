for sample in father mother proband; do
  samtools sort ${sample}_filtered.bam -o ${sample}_filtered_sorted.bam
  samtools index ${sample}_filtered_sorted.bam
done
