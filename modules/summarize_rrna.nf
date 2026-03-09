/*
========================================================================================
    SUMMARIZE_RRNA Module
========================================================================================
    Parses SortMeRNA and RiboDetector logs and generates a MultiQC generalstats
    custom content file with proper sample name mapping (FLI → sample_name).

    This avoids the MultiQC built-in SortMeRNA module, which derives sample names
    from reads file paths listed in the log (giving per-read-file rows instead of
    per-sample rows).
*/

process SUMMARIZE_RRNA {
    label 'process_single'
    executor 'local'

    input:
    path(sortmerna_logs)    // *.sortmerna.log files (or NO_SORTMERNA placeholder)
    path(ribodetector_logs) // *.ribodetector.log files (or NO_RIBODETECTOR placeholder)
    val(sample_info)        // List of maps: fli, sample_name, species

    output:
    path("rrna_pct_mqc.txt"), emit: summary

    script:
    def has_sortmerna    = sortmerna_logs.name    != 'NO_SORTMERNA'    ? 'True' : 'False'
    def has_ribodetector = ribodetector_logs.name != 'NO_RIBODETECTOR' ? 'True' : 'False'
    def sample_info_json = groovy.json.JsonOutput.toJson(sample_info)
    """
    #!/usr/bin/env python3
    import re, json, glob

    sample_info = json.loads('${sample_info_json}')
    fli_to_display = {
        s.get('fli', s['sample_name']): (
            s['sample_name'] + (' ({})'.format(s['species']) if s.get('species') else '')
        )
        for s in sample_info
    }

    has_sortmerna    = ${has_sortmerna}
    has_ribodetector = ${has_ribodetector}
    smr_data  = {}
    ribo_data = {}

    if has_sortmerna:
        for log_file in glob.glob("*.sortmerna.log"):
            fli     = log_file.replace(".sortmerna.log", "")
            display = fli_to_display.get(fli, fli)
            with open(log_file) as f:
                content = f.read()
            m = re.search(r'Total reads passing E-value threshold\\s*=\\s*\\d+\\s*\\(([\\d.]+)%?\\)', content)
            if m:
                smr_data[display] = float(m.group(1))

    if has_ribodetector:
        ansi_re = re.compile(r'\\x1b\\[[0-9;]*[A-Za-z]')
        for log_file in glob.glob("*.ribodetector.log"):
            fli     = log_file.replace(".ribodetector.log", "")
            display = fli_to_display.get(fli, fli)
            with open(log_file) as f:
                content = ansi_re.sub('', f.read())
            m_total = re.search(r'Processed\\s+(\\d+)\\s+sequences in total', content)
            m_rrna  = re.search(r'Detected\\s+(\\d+)\\s+rRNA sequences', content)
            if m_total and m_rrna:
                total = int(m_total.group(1))
                rrna  = int(m_rrna.group(1))
                if total > 0:
                    ribo_data[display] = rrna / total * 100

    with open("rrna_pct_mqc.txt", 'w') as f:
        f.write("# plot_type: 'generalstats'\\n")
        f.write("# id: 'rrna_pct'\\n")
        f.write("# pconfig:\\n")
        if has_sortmerna:
            f.write("#     pct_rrna_sortmerna:\\n")
            f.write("#         title: '% rRNA'\\n")
            f.write("#         description: 'Percent rRNA reads detected by SortMeRNA (subsampled reads)'\\n")
            f.write("#         namespace: 'SortMeRNA'\\n")
            f.write("#         suffix: '%'\\n")
            f.write("#         format: '{:,.1f}'\\n")
            f.write("#         min: 0\\n")
            f.write("#         max: 100\\n")
            f.write("#         scale: 'YlOrRd'\\n")
            f.write("#         placement: 95\\n")
        if has_ribodetector:
            f.write("#     pct_rrna_ribodetector:\\n")
            f.write("#         title: '% rRNA (RiboD)'\\n")
            f.write("#         description: 'Percent rRNA reads detected by RiboDetector (subsampled reads)'\\n")
            f.write("#         namespace: 'RiboDetector'\\n")
            f.write("#         suffix: '%'\\n")
            f.write("#         format: '{:,.1f}'\\n")
            f.write("#         min: 0\\n")
            f.write("#         max: 100\\n")
            f.write("#         scale: 'YlOrRd'\\n")
            f.write("#         placement: 96\\n")

        cols = ['Sample']
        if has_sortmerna:    cols.append('pct_rrna_sortmerna')
        if has_ribodetector: cols.append('pct_rrna_ribodetector')
        f.write('\\t'.join(cols) + '\\n')

        all_samples = sorted(set(list(smr_data.keys()) + list(ribo_data.keys())))
        for sample in all_samples:
            row = [sample]
            if has_sortmerna:
                row.append(f"{smr_data[sample]:.1f}" if sample in smr_data else '')
            if has_ribodetector:
                row.append(f"{ribo_data[sample]:.1f}" if sample in ribo_data else '')
            f.write('\\t'.join(row) + '\\n')

    print(f"rRNA summary: {len(smr_data)} SortMeRNA, {len(ribo_data)} RiboDetector samples")
    """
}
