/*
========================================================================================
    SUMMARIZE_RESULTS Module
========================================================================================
    Generates a consolidated summary table with key metrics from all QC modules:
    - FastQC: total reads, % duplicates, % GC, avg read length
    - Kraken2: % mtDNA, top genus, % top genus, top species, % top species
    - Sex determination: inferred sex, confidence
    - SortMeRNA: % rRNA (on subsampled reads)
    - RiboDetector: % rRNA (on subsampled reads)
*/

process SUMMARIZE_RESULTS {
    label 'process_single'
    executor 'local'
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path(fastqc_data)          // FastQC zip files
    path(kraken2_summary)      // kraken2_top_species_mqc.txt
    path(sex_summary)          // sex_determination_mqc.txt (or NO_SEX)
    path(sortmerna_logs)       // *.sortmerna.log files (or NO_SORTMERNA)
    path(ribodetector_logs)    // *.ribodetector.log files (or NO_RIBODETECTOR)
    path(rrna_kraken2_reports) // *_rrna.kraken2.report.txt files (or NO_RRNA_KRAKEN2)
    val(sample_info)           // List of maps with: sample_name, species

    output:
    path("qc_summary.tsv"), emit: summary

    script:
    def has_sex           = sex_summary.name           != 'NO_SEX'           ? 'True' : 'False'
    def has_sortmerna     = sortmerna_logs.name         != 'NO_SORTMERNA'    ? 'True' : 'False'
    def has_ribodetector  = ribodetector_logs.name      != 'NO_RIBODETECTOR' ? 'True' : 'False'
    def has_rrna_kraken2  = rrna_kraken2_reports.name   != 'NO_RRNA_KRAKEN2' ? 'True' : 'False'
    def sample_info_json  = groovy.json.JsonOutput.toJson(sample_info)
    """
    #!/usr/bin/env python3
    import os
    import re
    import json
    import zipfile
    import glob
    from collections import defaultdict

    # Sample metadata
    sample_info = json.loads('${sample_info_json}')
    sample_species = {s['sample_name']: s.get('species', '') for s in sample_info}

    # Map from FLI (flowcell ID) to sample_name
    fli_to_sample = {s.get('fli', s['sample_name']): s['sample_name'] for s in sample_info}

    # Data storage - keyed by sample_name
    data = defaultdict(dict)

    # Parse FastQC data from zip files
    for zip_file in glob.glob("*.zip"):
        try:
            with zipfile.ZipFile(zip_file, 'r') as zf:
                # Find the fastqc_data.txt file inside
                for name in zf.namelist():
                    if name.endswith('fastqc_data.txt'):
                        with zf.open(name) as f:
                            content = f.read().decode('utf-8')

                            # Extract sample name from zip filename
                            # Format: FLI_1_fastqc.zip or FLI_2_fastqc.zip
                            sample_base = zip_file.replace('_fastqc.zip', '')
                            # Remove _1 or _2 suffix to get FLI
                            if sample_base.endswith('_1') or sample_base.endswith('_2'):
                                read_num = sample_base[-1]
                                fli = sample_base[:-2]
                            else:
                                fli = sample_base
                                read_num = '1'

                            # Map FLI to sample_name
                            sample_name = fli_to_sample.get(fli, fli)

                            # Parse basic stats
                            in_basic_stats = False
                            for line in content.split('\\n'):
                                if line.startswith('>>Basic Statistics'):
                                    in_basic_stats = True
                                elif line.startswith('>>END_MODULE'):
                                    in_basic_stats = False
                                elif in_basic_stats and '\\t' in line:
                                    parts = line.split('\\t')
                                    if len(parts) >= 2:
                                        key, value = parts[0], parts[1]
                                        if key == 'Total Sequences':
                                            # Sum reads from R1 and R2
                                            prev = data[sample_name].get('total_reads', 0)
                                            data[sample_name]['total_reads'] = prev + int(value)
                                        elif key == '%GC' and read_num == '1':
                                            data[sample_name]['percent_gc'] = value
                                        elif key == 'Sequence length' and read_num == '1':
                                            data[sample_name]['read_length'] = value

                            # Parse duplication levels
                            in_dup = False
                            for line in content.split('\\n'):
                                if line.startswith('>>Sequence Duplication Levels'):
                                    in_dup = True
                                elif line.startswith('>>END_MODULE'):
                                    in_dup = False
                                elif in_dup and line.startswith('#Total Deduplicated'):
                                    parts = line.split('\\t')
                                    if len(parts) >= 2 and read_num == '1':
                                        dedup_pct = float(parts[1])
                                        data[sample_name]['percent_duplicates'] = f"{100 - dedup_pct:.1f}"
        except Exception as e:
            print(f"Warning: Could not parse {zip_file}: {e}")

    # Parse Kraken2 summary
    if "${kraken2_summary}" != "NO_KRAKEN2" and os.path.exists("${kraken2_summary}"):
        with open("${kraken2_summary}", 'r') as f:
            header = None
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\\t')
                if header is None:
                    header = parts
                    continue
                if len(parts) >= 6:
                    # Sample name might have species in parentheses, extract base name
                    sample_full = parts[0]
                    sample_name = sample_full.split(' (')[0] if ' (' in sample_full else sample_full

                    data[sample_name]['percent_mtdna'] = parts[1]
                    data[sample_name]['top_genus'] = parts[2]
                    data[sample_name]['percent_top_genus'] = parts[3]
                    data[sample_name]['top_species'] = parts[4]
                    data[sample_name]['percent_top_species'] = parts[5]

    # Parse sex determination summary if available
    has_sex = ${has_sex}
    if has_sex and "${sex_summary}" != "NO_SEX" and os.path.exists("${sex_summary}"):
        with open("${sex_summary}", 'r') as f:
            header = None
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\\t')
                if header is None:
                    header = parts
                    continue
                if len(parts) >= 2:
                    sample_full = parts[0]
                    sample_name = sample_full.split(' (')[0] if ' (' in sample_full else sample_full

                    data[sample_name]['inferred_sex'] = parts[1] if len(parts) > 1 else ''
                    data[sample_name]['sex_confidence'] = parts[2] if len(parts) > 2 else ''

    # Parse SortMeRNA logs
    # Log line: "Total reads passing E-value threshold = N (X.XX%)"
    # One log per FLI; average across FLIs for the same sample_name.
    has_sortmerna = ${has_sortmerna}
    if has_sortmerna:
        sortmerna_vals = defaultdict(list)
        for log_file in glob.glob("*.sortmerna.log"):
            fli = log_file.replace(".sortmerna.log", "")
            sample_name = fli_to_sample.get(fli, fli)
            with open(log_file) as f:
                content = f.read()
            m = re.search(r'Total reads passing E-value threshold\\s*=\\s*\\d+\\s*\\(([\\d.]+)%\\)', content)
            if m:
                sortmerna_vals[sample_name].append(float(m.group(1)))
        for sample_name, vals in sortmerna_vals.items():
            data[sample_name]['sortmerna_pct_rrna'] = f"{sum(vals)/len(vals):.2f}"

    # Parse RiboDetector logs
    # Log line: "rRNA sequences: N (X.XX%)"
    has_ribodetector = ${has_ribodetector}
    if has_ribodetector:
        ribodetector_vals = defaultdict(list)
        for log_file in glob.glob("*.ribodetector.log"):
            fli = log_file.replace(".ribodetector.log", "")
            sample_name = fli_to_sample.get(fli, fli)
            with open(log_file) as f:
                content = f.read()
            m = re.search(r'rRNA sequences:\\s*\\d+\\s*\\(([\\d.]+)%\\)', content, re.IGNORECASE)
            if m:
                ribodetector_vals[sample_name].append(float(m.group(1)))
        for sample_name, vals in ribodetector_vals.items():
            data[sample_name]['ribodetector_pct_rrna'] = f"{sum(vals)/len(vals):.2f}"

    # Parse rRNA Kraken2 reports
    # Report files are named {fli}_rrna.kraken2.report.txt
    has_rrna_kraken2 = ${has_rrna_kraken2}
    if has_rrna_kraken2:
        for report_file in glob.glob("*_rrna.kraken2.report.txt"):
            fli = report_file.replace("_rrna.kraken2.report.txt", "")
            sample_name = fli_to_sample.get(fli, fli)
            unclassified_pct = 0
            species_rows = []
            with open(report_file) as f:
                for line in f:
                    parts = line.strip().split('\\t')
                    if len(parts) < 6:
                        continue
                    rank = parts[3].strip()
                    pct = float(parts[0].strip())
                    name = parts[5].strip()
                    if rank == 'U':
                        unclassified_pct = pct
                    elif rank == 'S':
                        species_rows.append((pct, name))
            if species_rows:
                classified_pct = 100 - unclassified_pct
                top_pct, top_name = max(species_rows)
                pct_of_classified = (top_pct / classified_pct * 100) if classified_pct > 0 else 0
                data[sample_name]['rrna_top_species'] = top_name
                data[sample_name]['rrna_pct_top_species'] = f"{pct_of_classified:.1f}"

    # Write output TSV
    columns = [
        'sample_name', 'expected_species', 'total_reads', 'read_length',
        'percent_gc', 'percent_duplicates',
        'percent_mtdna', 'top_genus', 'percent_top_genus',
        'top_species', 'percent_top_species'
    ]
    if has_sex:
        columns.extend(['inferred_sex', 'sex_confidence'])
    if has_sortmerna:
        columns.append('sortmerna_pct_rrna')
    if has_ribodetector:
        columns.append('ribodetector_pct_rrna')
    if has_rrna_kraken2:
        columns.extend(['rrna_top_species', 'rrna_pct_top_species'])

    with open('qc_summary.tsv', 'w') as f:
        f.write('\\t'.join(columns) + '\\n')

        for sample_name in sorted(data.keys()):
            row = [
                sample_name,
                sample_species.get(sample_name, ''),
                str(data[sample_name].get('total_reads', '')),
                data[sample_name].get('read_length', ''),
                data[sample_name].get('percent_gc', ''),
                data[sample_name].get('percent_duplicates', ''),
                data[sample_name].get('percent_mtdna', ''),
                data[sample_name].get('top_genus', ''),
                data[sample_name].get('percent_top_genus', ''),
                data[sample_name].get('top_species', ''),
                data[sample_name].get('percent_top_species', '')
            ]
            if has_sex:
                row.extend([
                    data[sample_name].get('inferred_sex', ''),
                    data[sample_name].get('sex_confidence', '')
                ])
            if has_sortmerna:
                row.append(data[sample_name].get('sortmerna_pct_rrna', ''))
            if has_ribodetector:
                row.append(data[sample_name].get('ribodetector_pct_rrna', ''))
            if has_rrna_kraken2:
                row.extend([
                    data[sample_name].get('rrna_top_species', ''),
                    data[sample_name].get('rrna_pct_top_species', '')
                ])
            f.write('\\t'.join(row) + '\\n')

    print(f"Generated summary for {len(data)} samples")
    """
}
