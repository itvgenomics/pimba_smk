rule prepare_tax_assignments:
    input:
        tax_assignments = os.path.join(output_dir, 'output', file_name_raw + '_otus_tax_assignments.txt') if strategy == "otu" else os.path.join(output_dir, 'output', file_name_raw + '_asvs_tax_assignments.txt')
    output:
        prepared_tax_assignments = os.path.join(output_dir, 'output', 'plot_' + file_name_raw + '_tax_assignments.txt')
    shell:
        """
        cp {input.tax_assignments} {output.prepared_tax_assignments}
        
        COLUMNS=$(head -1 {output.prepared_tax_assignments} | sed -e 's/[^\t]//g' | wc -c)
        temp_file=$(mktemp)    
        if [ "$COLUMNS" -gt 3 ]; 
        then 
            echo -e "otu_id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tsimilarity\taux" > "$temp_file"
        else
            echo -e "otu_id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tsimilarity" > "$temp_file"
        fi

        cat "{output.prepared_tax_assignments}" >> "$temp_file"
        mv "$temp_file" "{output.prepared_tax_assignments}"

        sed -i 's/;/\t/g' {output.prepared_tax_assignments}
        sed -i 's/k__//g' {output.prepared_tax_assignments}
        sed -i 's/p__//g' {output.prepared_tax_assignments}
        sed -i 's/c__//g' {output.prepared_tax_assignments}
        sed -i 's/o__//g' {output.prepared_tax_assignments}
        sed -i 's/f__//g' {output.prepared_tax_assignments}
        sed -i 's/g__//g' {output.prepared_tax_assignments}
        sed -i 's/s__//g' {output.prepared_tax_assignments}
        sed -i 's/D_0__//g' {output.prepared_tax_assignments}
        sed -i 's/D_1__//g' {output.prepared_tax_assignments}
        sed -i 's/D_2__//g' {output.prepared_tax_assignments}
        sed -i 's/D_3__//g' {output.prepared_tax_assignments}
        sed -i 's/D_4__//g' {output.prepared_tax_assignments}
        sed -i 's/D_5__//g' {output.prepared_tax_assignments}
        sed -i 's/D_6__//g' {output.prepared_tax_assignments}
        sed -i 's/Unassigned/Unassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned\tUnassigned/g' {output.prepared_tax_assignments}
        """