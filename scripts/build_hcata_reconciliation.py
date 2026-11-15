#!/usr/bin/env python3
"""Item #17: Reconcile all 76 HCATA (PMC12511636, Bartz et al. 2025) Table S1 and 11 Table S2 studies."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import cellxgene_census
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.followup_organs import write_clean_csv

PMIDS_JSON = REPO.parent / "agescape-revision-data" / "tmp_organs" / "meta_probes" / "hcata_pmids.json"

# Exact rows from Table S1 of 42003_2025_8845_MOESM1_ESM.pdf (pages 2-4)
TABLE_S1 = [
    ("28965763", "Enge", "2017", "human", "0-54", "both", "Smart-Seq2", "No"),
    ("30093597", "Young", "2018", "human", "1-72", "unknown", "10x 3' v2", "No"),
    ("30315278", "Guo", "2018", "human", "1-25", "male", "10x 3' v2", "No"),
    ("30348985", "MacParland", "2018", "human", "21-66", "both", "10x 3' v2", "No"),
    ("30449713", "Wu", "2018", "human", "62-62", "male", "10x 3' v2", "No"),
    ("30518681", "Oetjen", "2018", "human", "24-84", "both", "10x 3' v2", "Yes"),
    ("30554520", "Reyfman", "2019", "human", "21-72", "both", "10x 3' v2", "No"),
    ("30566875", "Henry", "2018", "human", "25-31", "male", "10x 3' v2", "No"),
    ("30865899", "Fang", "2019", "human", "21-61", "both", "Makosco Drop-Seq", "No"),
    ("31023895", "Merrick", "2019", "human", "31-31", "female", "10x 3' v2", "No"),
    ("31178122", "Welch", "2019", "human", "18-75", "both", "10x 3' v2", "Yes"),
    ("31316211", "Schirmer", "2019", "human", "34-82", "both", "10x 3' v2", "Yes"),
    ("31436334", "Lukowski", "2019", "human", "42-80", "both", "10x 3' v2", "No"),
    ("31506348", "Wilson", "2019", "human", "52-74", "both", "10x 5' v2", "No"),
    ("31896769", "Liao", "2020", "human", "57-65", "both", "10x 3' v2", "No"),
    ("31924475", "Cillo", "2020", "human", "28-80", "both", "10x 3' v2", "Yes"),
    ("31980621", "Durante", "2020", "human", "44-86", "both", "10x 5' v2", "No"),
    ("32079746", "Park", "2020", "human", "0-40", "both", "10x 3' v2", "No"),
    ("32214235", "Han", "2020", "human", "0-83", "both", "microwell-seq", "No"),
    ("32327715", "Sole-Boldo", "2020", "human", "25-70", "male", "10x 3' v2", "No"),
    ("32451460", "Lee", "2020", "human", "50-86", "both", "10x 3' v2", "No"),
    ("32726565", "Deprez", "2020", "human", "24-35", "both", "10x 3' v2", "No"),
    ("32832599", "Adams", "2020", "human", "20-80", "both", "10x 3' v2", "Yes"),
    ("32888429", "Drokhlyansky", "2020", "human", "35-90", "both", "10x 3' v2", "Yes"),
    ("32948672", "Pappalardo", "2020", "human", "38-60", "both", "10x 3' v3", "No"),
    ("32971526", "Litvinukova", "2020", "human", "50-75", "both", "10x 3' v3", "No"),
    ("33017217", "Li", "2020", "human", "56-78", "both", "10x 3' v3", "No"),
    ("33192286", "Alsema", "2020", "human", "72-103", "both", "Modified Smart-Seq2", "No"),
    ("33212097", "Tosti", "2021", "human", "1.0-77", "both", "10x 3' v3 or v2", "No"),
    ("33287869", "He", "2020", "human", "31-31", "male", "10x 5'", "No"),
    ("33361824", "Trump", "2021", "human", "24-91", "both", "10x 3' v3", "No"),
    ("33397933", "Li", "2021", "human", "51-86", "female", "10x 3' v2", "No"),
    ("33558695", "Maynard", "2021", "human", "30-46", "both", "10x 3' v3", "No"),
    ("33579751", "King", "2021", "human", "3-9", "both", "10x 5' v2", "No"),
    ("33682422", "Sim", "2021", "human", "4-42", "both", "10x 3' v3", "No"),
    ("33861994", "Krishna", "2021", "human", "50-69", "both", "10x 5' v2", "No"),
    ("33915569", "Delorey", "2021", "human", "30-89", "both", "10x 3' v3", "No"),
    ("34025588", "Yang", "2021", "human", "28-39", "female", "10x 3' v3", "No"),
    ("34038734", "Olalekan", "2021", "human", "46-71", "both", "Makosco Drop-Seq", "No"),
    ("34051145", "Ayhan", "2021", "human", "24-60", "both", "10x 3' v2", "No"),
    ("34129837", "Williams", "2021", "human", "21-61", "both", "10x 3' v3", "No"),
    ("34400611", "Gan", "2021", "human", "18-31", "male", "10x 3' v2", "No"),
    ("34428183", "Shrestha", "2021", "human", "14-66", "both", "10x 3' v2", "Yes"),
    ("34447367", "Xu", "2021", "human", "60-77", "both", "10x 5'", "No"),
    ("34450029", "Pelka", "2021", "human", "35-91", "both", "10x 3' v2", "Yes"),
    ("34491790", "Yasen", "2021", "human", "27-56", "male", "10x 5'", "No"),
    ("34497389", "Elmentaite", "2021", "human", "1-29", "both", "10x 5' v2", "No"),
    ("34601654", "Rao", "2021", "human", "40-62", "male", "10x 3'", "No"),
    ("34702715", "Lee", "2022", "human", "9-55", "both", "10x", "No"),
    ("34722263", "Su", "2021", "human", "33-50", "both", "10x 3' v3", "No"),
    ("34990409", "Oliva", "2022", "human", "51-78", "both", "10x 3' v3", "No"),
    ("35021063", "Guilliams", "2022", "human", "29-75", "both", "10x 3' v3", "Yes"),
    ("35046994", "Tonglin", "2022", "human", "21-56", "both", "10x 5' v2", "No"),
    ("35073399", "Mende", "2022", "human", "19-35", "both", "Smart-Seq2", "No"),
    ("35216673", "COMBAT Consortium", "2022", "human", "2.0-9", "both", "10x 5' v1", "No"),
    ("35320732", "Ulrich", "2022", "human", "30-52", "female", "10x 3' v3", "No"),
    ("35389779", "Yazar", "2022", "human", "19-97", "both", "10x 3' v2", "No"),
    ("35389781", "Perez", "2022", "human", "20-83", "both", "10x 3' v2", "No"),
    ("35483358", "Patel", "2022", "human", "0-19", "both", "10x 3' v3", "No"),
    ("35538548", "Khaliq", "2022", "human", "37-83", "both", "10x 5'", "No"),
    ("35549404", "TS Consortium", "2022", "human", "22-74", "both", "Smart-Seq2, 10x 3' v3", "No"),
    ("35617956", "Gray", "2022", "human", "25-65", "female", "10x 3' v2", "No"),
    ("35688160", "Verhoeven", "2022", "human", "0-11.5", "unknown", "10x 3' v2", "No"),
    ("35858321", "van Zyl", "2022", "human", "29-75", "both", "10x 3' v3", "No"),
    ("35914526", "Fink", "2022", "human", "65-85", "both", "10x 3' v3", "No"),
    ("35959412", "Koenig", "2022", "human", "8.0-75.0", "both", "10x 5' v1", "No"),
    ("36108172", "Watanabe", "2022", "human", "33-81", "both", "10x 3' v2", "No"),
    ("36318921", "Herring", "2022", "human", "0-40", "both", "10x 3' v3", "No"),
    ("36543915", "Madissoon", "2023", "human", "29-75", "both", "10x 3' v3", "Yes"),
    ("36516485", "Perez", "2022", "human", "19-90", "both", "10x 5' v1", "Yes"),
    ("38198537", "Wang", "2024", "human", "69-93", "both", "10x 3' v3", "No"),
    ("26691212", "Li", "2016", "human", "37-37", "male", "Smart-Seq2", "No"),
    ("39578560", "Jin", "2025", "human", "23-54", "female", "10x 3' v3", "No"),
    ("https://explore.data.humancellatlas.org/projects/cc95ff89-2e68-4a08-a234-480eca21ce79", "Li", "2022", "human", "0-52", "both", "10x 3' v2", "No"),
    ("https://doi.org/10.1101/2022.02.09.479667", "Wang", "NA", "human", "53-65", "both", "10x 3' v3", "No"),
    ("https://doi.org/10.1101/2022.02.14.480397", "Zhu", "NA", "human", "63-96", "both", "10x 3' v3", "No"),
]

# Dispositions & primary-source reasons for all 76 rows
DISPOSITIONS = {
    "28965763": ("already_done_agescape", "Organ 04 Pancreas Enge (GSE81547) already analyzed in agescape followup"),
    "30093597": (
        "ineligible_sparse_and_disease",
        "PMC6104812 ('macroscopically normal regions of kidneys resected due to cancer (n=10) or for transplantation (n=2)'); ALL 40,268 Census cells have is_primary_data=False overlapping verbatim with Stewart 2019 (PMC7343525, Census 120e86b4); only 4 adult non-cancer transplant TxK1-4 (44F, 53F, 64M, 72M) + 5 RCC/VHL/pRCC tumor-resection margins + pediatric (<6 clean adult non-cancer donors)",
    ),
    "30315278": ("ineligible_sparse_donors", "Adult testis 1-25 yr: only 3 adult donors (17-25 yr)"),
    "30348985": ("ineligible_sparse_donors", "Normal liver GSE115469: n=5 donors (21-66 yr, <6 minimum)"),
    "30449713": ("ineligible_single_age", "Human kidney organoids/single 62-year-old donor (no age span)"),
    "30518681": ("blocked_missing_metadata", "GSE120221 bone marrow (n=20 healthy adults 24-84y in Oetjen Table S1/DS1): live GEO check confirms Matrix Market counts only; no public per-cell annotations or sex metadata"),
    "30554520": ("overlaps_hlca_and_already_done", "Lung IPF/controls (GSE122960): normal donors already integrated in HLCA (Lung_Integrated) already analyzed in agescape"),
    "30566875": ("ineligible_sparse_donors", "Normal prostate (Census): n=3 donors aged 25-31 yr (<6 minimum)"),
    "30865899": ("ineligible_disease_sparse", "Pancreas GSE101207 (T2D + normal controls): n=9 mixed T2D/healthy (<6 clean normal across adult span)"),
    "31023895": ("ineligible_single_age", "Adipose single 31F donor (no age span)"),
    "31178122": ("excluded_wrong_species_in_hcata", "Table S2 misattribution: Welch et al. 2019 (LIGER method paper; GEO GSE126836) is 15 mouse BNST brain Drop-seq datasets, not 38 human brain donors"),
    "31316211": ("ineligible_controlled_and_disease", "Multiple sclerosis brain lesions vs controls (n=9 mixed MS/control, EGA controlled)"),
    "31436334": ("ineligible_sparse_donors", "Retina (Census): only 2 normal donors with adult exact ages"),
    "31506348": ("ineligible_sparse_donors", "Kidney diabetic nephropathy vs controls: 3 controls + 3 DKD (<6 controls)"),
    "31896769": ("ineligible_sparse_donors", "Kidney GSE131685: n=3 donors aged 57-65 yr (<6 minimum)"),
    "31924475": ("ineligible_disease_sparse", "HNSCC cancer PBMC/TIL (GSE139324): Table S2 n=32 pools 26 HNSCC cancer patients with only 6 healthy PBMC donors"),
    "31980621": ("ineligible_disease", "Uveal melanoma tumor atlas (disease cohort, no normal aging arm)"),
    "32079746": ("ineligible_pediatric_fetal", "Thymic development atlas (0-40 yr): 0 adult exact-age donors >=20 yr in Census"),
    "32214235": ("ineligible_sparse_exact_donors", "Han 2020 HCL (microwell-seq, GSE134355): Census adult (>=20y) exact-age donors max 4 in blood and <=3 per solid tissue"),
    "32327715": ("already_done_agescape", "Organ 01 Skin Sole-Boldo (GSE130973) already analyzed in agescape followup"),
    "32451460": ("ineligible_disease", "Colorectal cancer GSE132465 (tumor vs matched adjacent mucosa, no healthy donor cohort)"),
    "32726565": ("overlaps_hlca", "Healthy airways Deprez 2020: donors (D322-D372) already integrated in HLCA (Lung_Integrated) already analyzed in agescape"),
    "32832599": ("overlaps_hlca", "IPF/control lung GSE136831 (Adams 2020): overlaps HLCA Lung_Integrated (donor 1372C) + IPF disease mix"),
    "32888429": ("ineligible_cancer_adjacent_sparse", "Broad SCP1038 / Drokhlyansky 2020: 10 CRC cancer-adjacent normal colon resections (35-90y), not 59 healthy human donors as listed in Table S2"),
    "32948672": ("ineligible_disease", "CSF/PBMC T cells in multiple sclerosis vs controls (small disease cohort)"),
    "32971526": ("ineligible_binned_ages", "Litvinukova 2020 adult human heart: 5-year age brackets only (no exact continuous ages); already excluded"),
    "33017217": ("ineligible_disease_sparse", "GSE155468 ascending thoracic aortic aneurysm (8 TAA patients + 3 controls)"),
    "33192286": ("ineligible_disease_narrow_elderly", "Alsema 2020 microglia (72-103y): AD/elderly autopsy focus, small sorted microglia cohort"),
    "33212097": ("ineligible_controlled_neonatal", "Pancreas EGA controlled access + adult/neonatal mix"),
    "33287869": ("ineligible_single_donor", "He 2020 15-organ atlas: single 31M donor (no age span)"),
    "33361824": ("ineligible_disease_controlled", "Trump 2021 COVID-19 airway (EGAS00001004689 controlled + acute viral infection)"),
    "33397933": ("ineligible_disease", "Vaginal wall from severe pelvic organ prolapse surgery + controls"),
    "33558695": ("ineligible_sparse_donors", "DLPFC spatial transcriptomics Visium (Maynard 2021): n=3 donors"),
    "33579751": ("ineligible_pediatric", "Pediatric tonsil B-cell maturation (ages 3-9 yr)"),
    "33682422": ("ineligible_sparse_pediatric", "Heart maturation (Sim 2021, 4-42y): Census has only 3 adult exact-age donors (35-42y)"),
    "33861994": ("ineligible_disease", "Clear cell renal cell carcinoma (ccRCC) tumor immune atlas"),
    "33915569": ("ineligible_disease", "COVID-19 autopsy multi-tissue atlas (acute fatal infection)"),
    "34025588": ("ineligible_pregnancy_disease", "Gestational diabetes placenta atlas (28-39y female pregnancy tissue)"),
    "34038734": ("ineligible_disease", "Metastatic ovarian cancer omental tumor microenvironment"),
    "34051145": ("ineligible_disease", "GSE160189 epilepsy surgical hippocampus resections (0 normal donors in Census)"),
    "34129837": ("blocked_missing_ages", "Oral mucosa GSE164241 (periodontitis + healthy controls): no donor ages in GEO sample records"),
    "34400611": ("ineligible_sparse_donors", "Intervertebral disc GSE160756 (degenerative disc + <=4 donors 18-31y)"),
    "34428183": ("ineligible_sparse_donors", "Pancreas GSE183568 (Shrestha 2021): n=5 non-diabetic donors including pediatric (14-66y, <6 adults)"),
    "34447367": ("ineligible_sparse_disease", "Alzheimer disease PBMC GSE181279: 3 AD + 2 control donors (<6)"),
    "34450029": ("ineligible_cancer_adjacent", "Pelka 2021 CRC (GSE178341): 36 'normal' 35-91y samples are tumor-adjacent normal colon mucosa from colorectal cancer resections"),
    "34491790": ("ineligible_disease", "Echinococcosis hepatic lesions (parasitic liver infection)"),
    "34497389": ("ineligible_fetal_pediatric", "Intestinal tract 1-29y (Elmentaite 2021): 0 adult exact-age donors in Census"),
    "34601654": ("ineligible_disease", "Heart failure GSE145154 (dilated/ischemic cardiomyopathy patients vs limited donor hearts)"),
    "34702715": ("ineligible_disease", "Chronic pancreatitis immune atlas (9-55y disease surgery)"),
    "34722263": ("ineligible_disease", "Renal cell carcinoma tumor atlas"),
    "34990409": ("ineligible_sparse_donors", "Olfactory mucosa GSE184117 (presbyosmia patients + 3 normosmia controls)"),
    "35021063": ("eligible_sorted_execute", "Guilliams 2022 Liver (GSE192742 / CZ e84f2780): 18 normal exact-age donors overall in Census (10 at 10x 3' v3 cell, 28-77y; FACS CD45+/CD45- sorted -> gene + P1 only)"),
    "35046994": ("ineligible_disease", "Aplastic anemia bone marrow scRNA-seq"),
    "35073399": ("ineligible_young_sparse", "Healthy adult HSPC (19-35y, Smart-Seq2): narrow young span + sorted HSPCs"),
    "35216673": ("ineligible_missing_ages_and_disease", "COMBAT Consortium COVID-19 blood atlas: 0 exact donor ages in Census + acute sepsis/COVID"),
    "35320732": ("ineligible_sparse_disease", "Fallopian tube GSE178101 (4 normal donors 30-52y + hydrosalpinx)"),
    "35389779": ("eligible_open_execute", "Yazar 2022 OneK1K PBMC (GSE196830 / CZ 3faad104): 981 normal donors (976 aged >=20, 19-97y, 565F/416M), open exact age, NOT selected by HCATA despite top power"),
    "35389781": ("eligible_open_execute", "Perez 2022 SLE/Control PBMC (GSE174188 / CZ 218acb0f): 99 normal control donors (21-75y exact age, 97F/2M; SLE excluded), NOT selected by HCATA"),
    "35483358": ("ineligible_pediatric_disease", "Rhabdomyosarcoma pediatric atlas (ages 0-19y)"),
    "35538548": ("ineligible_disease", "Colorectal cancer single-cell stratification"),
    "35549404": ("already_covered_ts2", "Tabula Sapiens Consortium (TS1/TS2): TS2 collection e5f58829 audited and analyzed in item #15"),
    "35617956": ("ineligible_sparse_and_brca", "Breast atlas GSE180878: mixed BRCA mutation carriers with only 3 WT non-carrier controls >=6 threshold not met"),
    "35688160": ("ineligible_pediatric_disease", "Neuroblastoma pediatric immune atlas (0-11.5y)"),
    "35858321": ("ineligible_sparse_donors", "Ocular anterior segment GSE199013: 6 adult exact-age donors overall in Census, max 5 per individual structure (<6)"),
    "35914526": ("ineligible_cancer_adjacent_elderly", "Human ureter GSE184111: 10 normal margins from bladder-cancer radical cystectomies, clustered 65-85y (narrow surgical cancer cohort)"),
    "35959412": ("overlaps_gao_heart", "Koenig 2022 donor vs DCM heart: already incorporated inside Item #14 Gao Heart GSE290367 and already excluded/separated per #14 rules"),
    "36108172": ("ineligible_cancer_adjacent_sparse", "COPD/normal lung GSE173896 (Watanabe 2022): 7 normal donors from primary lung cancer lobectomy resections (4 smokers + 3 never-smokers) -> surgical cancer margins"),
    "36318921": ("ineligible_developmental_sparse", "Herring 2022 PFC (GSE168408): 0-40y developmental series dominated by fetal/neonatal/adolescent brains"),
    "36543915": ("overlaps_hlca_mixed_platforms", "Madissoon 2023 lung/airway: overlaps HLCA (290B, 292B, 296C, 298C), mixes 5 airway regions and 3 platforms with n=7 exact donors; Table S1 pasted van Zyl eye cell types onto row"),
    "36516485": ("already_done_agescape", "Perez 2022 skeletal muscle (age 19-90): already analyzed in agescape muscle followup"),
    "38198537": ("ineligible_narrow_elderly_disease", "Wang 2024 substantia nigra GSE184950 (69-93y): glioma/elderly cohort"),
    "26691212": ("ineligible_single_donor", "Li 2016 pancreas: single 37M donor (no age span)"),
    "39578560": ("already_done_agescape", "Jin 2025 human ovary: already analyzed in agescape ovary followup"),
    "https://explore.data.humancellatlas.org/projects/cc95ff89-2e68-4a08-a234-480eca21ce79": ("ineligible_pooled_no_demux", "Li 2022 HCA Bone Marrow/PBMC: pooled genetically mixed libraries without per-cell donor age demultiplexing"),
    "https://doi.org/10.1101/2022.02.09.479667": ("ineligible_disease_transplant", "Wang bioRxiv 2022 liver transplant early allograft dysfunction (disease/transplant surgery, NA year)"),
    "https://doi.org/10.1101/2022.02.14.480397": ("ineligible_disease", "Zhu bioRxiv 2022 Parkinson disease brain (disease surgery/autopsy cohort, NA year)"),
}

# Exact Table S2 rows (page 5 of 42003_2025_8845_MOESM1_ESM.pdf)
TABLE_S2 = [
    ("30518681", "Oetjen 2018", 24, 84, 25, "Bone Marrow", "20 healthy adults (24-84y); live GEO check confirms Matrix Market files only without public per-cell annotations or sex", "blocked_missing_metadata"),
    ("31178122", "Welch 2019", 18, 66, 38, "Brain", "Misattribution: GSE126836 is 15 mouse BNST Drop-seq brain samples + LIGER paper, not 38 human brain donors", "excluded_wrong_species_in_hcata"),
    ("31316211", "Schirmer 2019", 34, 82, 9, "Brain", "Multiple sclerosis brain lesions vs controls (n=9 total mixed MS/control, EGA controlled)", "ineligible_controlled_and_disease"),
    ("31924475", "Cillo 2020", 29, 80, 32, "Blood", "GSE139324: pools 26 HNSCC head-and-neck cancer patients with only 6 healthy PBMC donors", "ineligible_disease_sparse"),
    ("32832599", "Adams 2020", 20, 67, 28, "Lung", "GSE136831: overlaps HLCA Lung_Integrated (donor 1372C already analyzed) + IPF disease mix", "overlaps_hlca"),
    ("32888429", "Drokhlyansky 2020", 35, 90, 59, "Colon", "Broad SCP1038: 10 CRC cancer-adjacent normal colon donors (35-90y), not 59 healthy human donors", "ineligible_cancer_adjacent_sparse"),
    ("34428183", "Shrestha 2021", 14, 66, 14, "Pancreas", "GSE183568: only 5 non-diabetic donors (including pediatric) mixed with T1D/AAB+ donors", "ineligible_sparse_donors"),
    ("34450029", "Pelka 2021", 35, 90, 36, "Colon", "GSE178341: 36 'normal' 35-91y samples are tumor-adjacent normal mucosa from CRC surgical resections", "ineligible_cancer_adjacent"),
    ("35021063", "Guilliams 2022", 28, 72, 7, "Liver", "GSE192742 / CZ e84f2780: Census contains 18 normal exact-age donors overall (10 at 10x 3' v3 cell, 28-77y; FACS sorted CD45+/CD45-)", "eligible_sorted_execute"),
    ("36516485", "Perez 2022", 19, 90, 17, "Skeletal Muscle", "Already analyzed in agescape muscle followup", "already_done_agescape"),
    ("36543915", "Madissoon 2023", 29, 64, 22, "Lung", "Overlaps HLCA (290B, 292B, 296C, 298C), mixes 5 airway regions and 3 platforms (n=7 exact donors); Table S1 erroneously pasted eye cell types", "overlaps_hlca_mixed_platforms"),
]


def load_pmid_metadata() -> dict[str, dict[str, str]]:
    if PMIDS_JSON.exists():
        hp = json.loads(PMIDS_JSON.read_text(encoding="utf-8"))
        return {str(x["pmid"]): x for x in hp}
    existing = REPO / "results" / "followup-hcata-pmc12511636-studies.csv"
    if existing.exists():
        df = pd.read_csv(existing, dtype=str).fillna("")
        return {
            row["identifier"]: {"pmc": row["pmc"], "doi": row["doi"], "title": row["title"]}
            for _, row in df.iterrows()
        }
    return {}


def main() -> None:
    pmid_meta = load_pmid_metadata()

    with cellxgene_census.open_soma(census_version="2025-11-08") as census:
        ds = census["census_info"]["datasets"].read().concat().to_pandas()
    ds["doi_low"] = ds["collection_doi"].fillna("").str.lower()
    doi_to_census = {}
    for doi, grp in ds[ds["doi_low"] != ""].groupby("doi_low"):
        doi_to_census[doi] = {
            "census_collection_id": str(grp["collection_id"].iloc[0]),
            "census_n_datasets": int(len(grp)),
            "census_total_cells": int(grp["dataset_total_cell_count"].sum()),
        }

    s1_rows = []
    for idx, (ident, author, year, species, age_range, sex, platform, selected) in enumerate(TABLE_S1, 1):
        meta = pmid_meta.get(ident, {})
        doi = str(meta.get("doi", "")).lower()
        cinfo = doi_to_census.get(doi, {})
        disp, reason = DISPOSITIONS[ident]
        s1_rows.append(
            {
                "table_s1_row": idx,
                "identifier": ident,
                "first_author": author,
                "year": year,
                "species": species,
                "hcata_age_range": age_range,
                "sex": sex,
                "platform": platform,
                "selected_for_deg_in_hcata": selected,
                "pmc": meta.get("pmc", ""),
                "doi": meta.get("doi", ""),
                "title": meta.get("title", ""),
                "census_collection_id": cinfo.get("census_collection_id", ""),
                "census_n_datasets": cinfo.get("census_n_datasets", 0),
                "census_total_cells": cinfo.get("census_total_cells", 0),
                "agescape_disposition": disp,
                "primary_source_verification_note": reason,
            }
        )

    s1_df = pd.DataFrame(s1_rows)
    assert len(s1_df) == 76, f"Expected 76 Table S1 rows, got {len(s1_df)}"
    write_clean_csv(s1_df, REPO / "results" / "followup-hcata-pmc12511636-studies.csv")

    s2_rows = []
    for idx, (pmid, author, ymin, ymax, n_samp, tissue, forensic, disp) in enumerate(TABLE_S2, 1):
        meta = pmid_meta.get(pmid, {})
        s2_rows.append(
            {
                "table_s2_row": idx,
                "pmid": pmid,
                "author_year": author,
                "doi": meta.get("doi", ""),
                "hcata_youngest_year": ymin,
                "hcata_oldest_year": ymax,
                "hcata_sample_size": n_samp,
                "hcata_tissue": tissue,
                "agescape_disposition": disp,
                "primary_source_forensic_note": forensic,
            }
        )
    s2_df = pd.DataFrame(s2_rows)
    assert len(s2_df) == 11, f"Expected 11 Table S2 rows, got {len(s2_df)}"
    write_clean_csv(s2_df, REPO / "results" / "followup-hcata-pmc12511636-table-s2-audit.csv")

    print("Disposition counts across 76 Table S1 studies:")
    print(s1_df["agescape_disposition"].value_counts().to_string())


if __name__ == "__main__":
    main()
