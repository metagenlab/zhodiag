process KRAKEN2BACTERIA_COMBINEKREPORTS {
  label 'process_single'

  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/krakentools:1.2--pyh5e36f6f_0':
      'biocontainers/krakentools:1.2--pyh5e36f6f_0' }"

  input:
  // NOTE: kreports is a LIST of paths (because of .collect()).
  // stageAs closure runs per file and we derive a unique staged name:
  path(kreports, stageAs: { p ->
      def db   = p.getParent().getParent().getFileName().toString()
      def conf = p.getParent().getFileName().toString()
      "${db}_${conf}_${p.getFileName().toString()}"
  })

  output:
  path("*/*/krakentools_combine_kreports.tsv"), emit: txt
  path("*/*/versions_krakentools.yml"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args    = task.ext.args ?: ''
  def VERSION = '1.2'

  // Use the first file to infer db/conf for the out_dir (same logic as before)
  def first_report = kreports[0]
  def fname = first_report.getFileName().toString()
  def m = fname =~ /^.+?_(.+?)_(conf\d+(?:\.\d+)?)_report\.txt$/
  if( !m.matches() ) throw new RuntimeException("Filename does not match expected format: ${fname}")

  def db   = m[0][1]
  def conf = m[0][2]
  def out_dir = "${db}/${conf}"

  """
  mkdir -p ${out_dir}

  combine_kreports.py \\
      -r ${kreports.collect{ it.getFileName().toString() }.join(' ')} \\
      -o ${out_dir}/krakentools_combine_kreports.tsv \\
      --display-headers \\
      ${args}

  cat <<-END_VERSIONS > ${out_dir}/versions_krakentools.yml
  "${task.process}":
      combine_kreports.py: ${VERSION}
  END_VERSIONS
  """
}
