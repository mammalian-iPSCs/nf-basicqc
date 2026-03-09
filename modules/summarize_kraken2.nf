/*
========================================================================================
    SUMMARIZE_KRAKEN2 Module
========================================================================================
    Parses Kraken2 reports and generates:
    - General Stats table with top species/genus per sample
    - Modified Kraken reports with unclassified removed and percentages recalculated
      (for interactive MultiQC Kraken plot without unclassified)
*/

process SUMMARIZE_KRAKEN2 {
    label 'process_single'
    executor 'local'

    input:
    path(reports)
    val(sample_info)  // List of maps with: sample_name, species

    output:
    path("kraken2_top_species_mqc.txt"), emit: summary
    path("*_classified.kraken2.report.txt"), emit: classified_reports

    script:
    // Build a lookup map from sample_name to display name (with species)
    def name_lookup = sample_info.collectEntries { info ->
        def display = info.sample_name
        if (info.species) {
            display = "${display} (${info.species})"
        }
        [(info.sample_name): display]
    }
    def name_lookup_json = groovy.json.JsonOutput.toJson(name_lookup)
    """
    #!/usr/bin/env python3
    import os
    import json
    import glob

    # Sample name lookup from samplesheet metadata
    name_lookup = json.loads('${name_lookup_json}')

    # Parse all kraken2 reports and extract data
    results = []

    for report_file in sorted(glob.glob("*.kraken2.report.txt")):
        # Extract sample ID from filename (remove .kraken2.report.txt)
        sample_id = report_file.replace(".kraken2.report.txt", "")

        # Get display name from lookup, fall back to sample_id
        display_name = name_lookup.get(sample_id, sample_id)

        top_species = "Unknown"
        top_species_percent = 0.0
        top_genus = "Unknown"
        top_genus_percent = 0.0
        percent_unclassified = 0.0

        # First pass: get unclassified percentage and find top taxa
        lines = []
        with open(report_file, 'r') as f:
            for line in f:
                lines.append(line)
                parts = line.strip().split('\\t')
                if len(parts) >= 6:
                    percent = float(parts[0].strip())
                    rank = parts[3].strip()
                    taxon = parts[5].strip()

                    if rank == 'U':
                        percent_unclassified = percent
                    elif rank == 'G' and percent > top_genus_percent:
                        top_genus_percent = percent
                        top_genus = taxon
                    elif rank == 'S' and percent > top_species_percent:
                        top_species_percent = percent
                        top_species = taxon

        # Calculate % classified
        percent_classified = 100.0 - percent_unclassified

        # Calculate top genus as % of classified reads
        top_genus_pct_of_classified = (top_genus_percent / percent_classified * 100.0) if percent_classified > 0 else 0

        results.append((display_name, sample_id, top_species, top_species_percent, top_genus, top_genus_pct_of_classified, percent_classified))

        # Create modified Kraken report without unclassified, with recalculated percentages
        # This will be used by MultiQC for the interactive plot
        output_file = f"{sample_id}_classified.kraken2.report.txt"
        with open(output_file, 'w') as f_out:
            for line in lines:
                parts = line.strip().split('\\t')
                if len(parts) >= 6:
                    rank = parts[3].strip()

                    # Skip unclassified line
                    if rank == 'U':
                        continue

                    # Recalculate percentage relative to classified reads
                    if percent_classified > 0:
                        orig_percent = float(parts[0].strip())
                        new_percent = (orig_percent / percent_classified) * 100.0
                        parts[0] = f"  {new_percent:.2f}"

                    f_out.write('\\t'.join(parts) + '\\n')

    # Write MultiQC custom content file as a table section
    with open("kraken2_top_species_mqc.txt", 'w') as f:
        f.write("# plot_type: 'table'\\n")
        f.write("# section_name: 'Kraken2 Species Identification'\\n")
        f.write("# description: 'Top species and genus identified from mtDNA database'\\n")
        f.write("# pconfig:\\n")
        f.write("#     id: 'kraken2_species_table'\\n")
        f.write("#     namespace: 'Kraken2'\\n")
        f.write("# headers:\\n")
        f.write("#     percent_classified:\\n")
        f.write("#         title: '% mtDNA'\\n")
        f.write("#         description: 'Percent of reads classified (mitochondrial)'\\n")
        f.write("#         suffix: '%'\\n")
        f.write("#         format: '{:,.2f}'\\n")
        f.write("#     top_genus:\\n")
        f.write("#         title: 'Top Genus'\\n")
        f.write("#         description: 'Most abundant genus detected'\\n")
        f.write("#     percent_top_genus:\\n")
        f.write("#         title: '% Top Genus'\\n")
        f.write("#         description: 'Percent of classified reads'\\n")
        f.write("#         suffix: '%'\\n")
        f.write("#         format: '{:,.1f}'\\n")
        f.write("#     top_species:\\n")
        f.write("#         title: 'Top Species'\\n")
        f.write("#         description: 'Most abundant species detected'\\n")
        f.write("#     percent_top_species:\\n")
        f.write("#         title: '% Top Species'\\n")
        f.write("#         description: 'Percent of total reads'\\n")
        f.write("#         suffix: '%'\\n")
        f.write("#         format: '{:,.2f}'\\n")
        f.write("Sample\\tpercent_classified\\ttop_genus\\tpercent_top_genus\\ttop_species\\tpercent_top_species\\n")

        for display_name, sample_id, top_species, top_species_percent, top_genus, top_genus_pct_classified, percent_classified in results:
            f.write(f"{display_name}\\t{percent_classified:.2f}\\t{top_genus}\\t{top_genus_pct_classified:.1f}\\t{top_species}\\t{top_species_percent:.2f}\\n")

    print(f"Processed {len(results)} Kraken2 reports")
    print(f"Created modified reports without unclassified reads")
    """
}

/*
========================================================================================
    SUMMARIZE_RRNA_KRAKEN2 Module
========================================================================================
    Parses rRNA Kraken2 (SILVA SSU) reports and generates a MultiQC custom content
    table with top species/genus per sample.

    Reports are named {fli}_rrna.kraken2.report.txt (one per FLI, not per sample_name).
    Results are aggregated by sample_name if multiple FLIs share a sample.
*/

process SUMMARIZE_RRNA_KRAKEN2 {
    label 'process_single'
    executor 'local'

    input:
    path(reports)
    val(sample_info)  // List of maps with: fli, sample_name, species

    output:
    path("rrna_kraken2_species_mqc.txt"), emit: summary

    script:
    def fli_to_display = sample_info.collectEntries { info ->
        def display = info.sample_name
        if (info.species) {
            display = "${display} (${info.species})"
        }
        [(info.get('fli', info.sample_name)): display]
    }
    def fli_to_display_json = groovy.json.JsonOutput.toJson(fli_to_display)
    """
    #!/usr/bin/env python3
    import json, glob
    from collections import defaultdict

    fli_to_display = json.loads('${fli_to_display_json}')

    # Aggregate results by display_name (in case multiple FLIs per sample)
    agg = defaultdict(lambda: {
        'classified_pct_sum': 0.0,
        'top_genus': None, 'top_genus_pct': 0.0,
        'top_species': None, 'top_species_pct': 0.0,
        'count': 0
    })

    for report_file in sorted(glob.glob("*_rrna.kraken2.report.txt")):
        fli = report_file.replace("_rrna.kraken2.report.txt", "")
        display = fli_to_display.get(fli, fli)

        unclassified_pct = 0.0
        genus_rows = []
        species_rows = []
        with open(report_file) as f:
            for line in f:
                parts = line.strip().split('\\t')
                if len(parts) < 6:
                    continue
                rank = parts[3].strip()
                pct  = float(parts[0].strip())
                name = parts[5].strip()
                if rank == 'U':
                    unclassified_pct = pct
                elif rank == 'G':
                    genus_rows.append((pct, name))
                elif rank == 'S':
                    species_rows.append((pct, name))

        classified_pct = 100.0 - unclassified_pct
        r = agg[display]
        r['classified_pct_sum'] += classified_pct
        r['count'] += 1

        if species_rows:
            top_pct, top_name = max(species_rows)
            pct_of_classified = (top_pct / classified_pct * 100) if classified_pct > 0 else 0
            if pct_of_classified > r['top_species_pct']:
                r['top_species'] = top_name
                r['top_species_pct'] = pct_of_classified

        if genus_rows:
            top_pct, top_name = max(genus_rows)
            pct_of_classified = (top_pct / classified_pct * 100) if classified_pct > 0 else 0
            if pct_of_classified > r['top_genus_pct']:
                r['top_genus'] = top_name
                r['top_genus_pct'] = pct_of_classified

    with open("rrna_kraken2_species_mqc.txt", 'w') as f:
        f.write("# plot_type: 'table'\\n")
        f.write("# section_name: 'rRNA Species ID (SILVA)'\\n")
        f.write("# description: 'Top species and genus from SILVA SSU Kraken2 on SortMeRNA rRNA reads'\\n")
        f.write("# pconfig:\\n")
        f.write("#     id: 'rrna_kraken2_species_table'\\n")
        f.write("#     namespace: 'rRNA Kraken2'\\n")
        f.write("# headers:\\n")
        f.write("#     rrna_pct_classified:\\n")
        f.write("#         title: '% Classified'\\n")
        f.write("#         description: 'Percent of rRNA reads classified against SILVA SSU'\\n")
        f.write("#         suffix: '%'\\n")
        f.write("#         format: '{:,.2f}'\\n")
        f.write("#     rrna_top_genus:\\n")
        f.write("#         title: 'Top Genus'\\n")
        f.write("#         description: 'Most abundant genus in rRNA reads (SILVA)'\\n")
        f.write("#     rrna_pct_top_genus:\\n")
        f.write("#         title: '% Top Genus'\\n")
        f.write("#         description: 'Percent of classified rRNA reads'\\n")
        f.write("#         suffix: '%'\\n")
        f.write("#         format: '{:,.1f}'\\n")
        f.write("#     rrna_top_species:\\n")
        f.write("#         title: 'Top Species'\\n")
        f.write("#         description: 'Most abundant species in rRNA reads (SILVA)'\\n")
        f.write("#     rrna_pct_top_species:\\n")
        f.write("#         title: '% Top Species'\\n")
        f.write("#         description: 'Percent of classified rRNA reads'\\n")
        f.write("#         suffix: '%'\\n")
        f.write("#         format: '{:,.2f}'\\n")
        f.write("Sample\\trrna_pct_classified\\trrna_top_genus\\trrna_pct_top_genus\\trrna_top_species\\trrna_pct_top_species\\n")

        for display, r in sorted(agg.items()):
            avg_classified = r['classified_pct_sum'] / r['count'] if r['count'] > 0 else 0
            f.write(
                f"{display}\\t{avg_classified:.2f}"
                f"\\t{r['top_genus'] or 'Unknown'}\\t{r['top_genus_pct']:.1f}"
                f"\\t{r['top_species'] or 'Unknown'}\\t{r['top_species_pct']:.2f}\\n"
            )

    print(f"rRNA Kraken2 summary: {len(agg)} samples")
    """
}
