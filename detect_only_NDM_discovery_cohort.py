import pandas as pd
from collections import defaultdict

def load_disease_mapping(cohort_file):
    """Load patient to disease mapping, focusing on NDM cases"""
    patient_disease = {}
    with open(cohort_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                patient = parts[0]
                category = parts[1]
                
                # Focus only on NDM cases (Solved and Unsolved)
                if "NDM" in category:
                    if category.startswith('Solved_NDM'):
                        patient_disease[patient] = 'Solved NDM'
                    else:
                        patient_disease[patient] = 'Unsolved NDM'
    return patient_disease

def process_shared_regions(shared_file, patient_disease):
    """Process shared regions focusing only on NDM patients"""
    regions_data = []
    
    with open(shared_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                chrom = parts[0]
                start = parts[1]
                end = parts[2]
                patients = [p for p in parts[5].split(',') if p in patient_disease]
                
                # Only process regions shared by NDM patients
                if len(patients) > 0:
                    # Get disease status for these patients
                    patient_diseases = []
                    disease_counts = defaultdict(int)
                    for p in patients:
                        disease = patient_disease[p]
                        patient_diseases.append(f"{p}({disease})")
                        disease_counts[disease] += 1
                    
                    regions_data.append({
                        'Chromosome': chrom,
                        'Start': start,
                        'End': end,
                        'Region': f"{chrom}:{start}-{end}",
                        'N_Patients': len(patients),
                        'Solved_NDM': disease_counts.get('Solved NDM', 0),
                        'Unsolved_NDM': disease_counts.get('Unsolved NDM', 0),
                        'Patients': ", ".join(patients),
                        'Patient_Diseases': ", ".join(patient_diseases),
                        'Status': 'Mixed' if len(disease_counts) > 1 else list(disease_counts.keys())[0]
                    })
    
    return pd.DataFrame(regions_data)

def save_ndm_reports(df):
    """Save NDM-specific reports"""
    # Main NDM shared regions
    df.to_csv('ndm_shared_regions.tsv', sep='\t', index=False)
    print("Main NDM report saved as ndm_shared_regions.tsv")
    
    # Regions shared between Solved and Unsolved cases
    mixed_regions = df[df['Status'] == 'Mixed']
    if len(mixed_regions) > 0:
        mixed_regions.to_csv('mixed_ndm_regions.tsv', sep='\t', index=False)
        print(f"Found {len(mixed_regions)} regions shared between Solved/Unsolved NDM cases")
        print("Saved as mixed_ndm_regions.tsv")
    else:
        print("No regions shared between Solved and Unsolved NDM cases found")

    # Summary counts
    summary = df.groupby('Status').agg({
        'Region': 'count',
        'N_Patients': 'sum'
    }).reset_index()
    summary.columns = ['Status', 'N_Regions', 'Total_Patients']
    summary.to_csv('ndm_summary.tsv', sep='\t', index=False)
    print("Summary statistics saved as ndm_summary.tsv")

if __name__ == "__main__":
    cohort_file = "new_cohort.txt"
    shared_regions_file = "unsolved_cohort_1000_8.tsv"
    
    try:
        print("Loading NDM patient information...")
        patient_disease = load_disease_mapping(cohort_file)
        print(f"Found {len(patient_disease)} NDM patients ({sum(1 for v in patient_disease.values() if v == 'Solved NDM')} Solved, "
              f"{sum(1 for v in patient_disease.values() if v == 'Unsolved NDM')} Unsolved)")
        
        print("\nProcessing shared regions...")
        ndm_df = process_shared_regions(shared_regions_file, patient_disease)
        
        if len(ndm_df) > 0:
            print(f"\nFound {len(ndm_df)} regions shared among NDM patients")
            save_ndm_reports(ndm_df)
            
            # Print key findings
            print("\nKEY FINDINGS")
            print("=" * 40)
            print(f"Total regions shared by NDM patients: {len(ndm_df)}")
            print(f"Regions exclusive to Solved NDM: {len(ndm_df[ndm_df['Status'] == 'Solved NDM'])}")
            print(f"Regions exclusive to Unsolved NDM: {len(ndm_df[ndm_df['Status'] == 'Unsolved NDM'])}")
            print(f"Regions shared between both: {len(ndm_df[ndm_df['Status'] == 'Mixed'])}")
            
            if len(ndm_df[ndm_df['Status'] == 'Mixed']) > 0:
                print("\nMost shared region between Solved/Unsolved NDM:")
                mixed = ndm_df[ndm_df['Status'] == 'Mixed'].sort_values('N_Patients', ascending=False).iloc[0]
                print(f"Region: {mixed['Region']}")
                print(f"Patients: {mixed['Patient_Diseases']}")
        else:
            print("No shared regions found among NDM patients")
            
    except FileNotFoundError as e:
        print(f"Error: File not found - {str(e)}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")