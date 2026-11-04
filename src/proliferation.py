import re


PROLIFERATIVE = "proliferative_capable"
POST_MITOTIC = "post_mitotic"
MEIOTIC_OR_MIXED = "meiotic_or_mixed"
UNCERTAIN = "uncertain"


def _has(label, pattern):
    return re.search(pattern, label) is not None


def classify_mammal(label):
    if _has(label, r"neuron|interneuron"):
        return POST_MITOTIC, "mature neuron"
    if label in {"cardiomyocyte", "skeletal muscle cell", "podocyte", "erythrocyte", "oligodendrocyte"}:
        return POST_MITOTIC, "terminally differentiated mature cell"
    if _has(label, r"spermatid|spermatozoa"):
        return POST_MITOTIC, "post-meiotic male germ cell"
    if _has(label, r"multiciliated|ciliated cell|ciliated columnar"):
        return POST_MITOTIC, "terminally differentiated ciliated cell"
    if label in {"pulmonary alveolar type 1 cell", "plasma cell", "neutrophil", "basophil", "granulocyte"}:
        return POST_MITOTIC, "terminally differentiated mature cell"
    if _has(label, r"spermatocyte"):
        return MEIOTIC_OR_MIXED, "meiotic germ-cell stage"
    if _has(label, r"stem|progenitor|precursor|erythroblast|proerythroblast|thymocyte|promonocyte|granulocytopoietic|plasmablast"):
        return PROLIFERATIVE, "stem, progenitor, or developing cell"
    if _has(label, r"\b[bt] cell\b|t cell|nk cell|nkt|natural killer|innate lymphoid|lymphocyte"):
        return PROLIFERATIVE, "lymphoid cell with activation-associated proliferation"
    if _has(label, r"monocyte|macrophage|microglial|dendritic cell|mast cell|langerhans cell"):
        return PROLIFERATIVE, "myeloid or innate immune cell with proliferative capacity"
    if _has(label, r"fibroblast|myofibroblast|endothelial|endocardial|fenestrated|pericyte|mesangial|stellate|smooth muscle|adventitial|leptomeningeal"):
        return PROLIFERATIVE, "stromal or vascular cell with proliferative capacity"
    if _has(label, r"hepatocyte|keratinocyte|urothelial|pneumocyte|acinar cell|pancreatic [abd] cell|mesothelial"):
        return PROLIFERATIVE, "epithelial or parenchymal cell with proliferative capacity"
    if _has(label, r"kidney|renal|proximal tubule|collecting duct|connecting tubule|distal convoluted|loop of henle|parietal epithelial"):
        return PROLIFERATIVE, "renal epithelial or vascular cell with proliferative capacity"
    if _has(label, r"basal cell|respiratory tract hillock|respiratory tract suprabasal|club cell|secretory cell|serous secreting|astrocyte|schwann cell|satellite cell|mesenchymal"):
        return PROLIFERATIVE, "facultative progenitor or proliferative-capable mature cell"
    if label in {"luminal epithelial cell of mammary gland", "kidney cortex artery cell"}:
        return PROLIFERATIVE, "epithelial or vascular cell with proliferative capacity"
    return UNCERTAIN, "mixed, generic, or insufficiently resolved annotation"


def classify_fly(label):
    if _has(label, r"intestinal stem cell|follicle stem cell|enteroblast|16-cell germline cyst|young germ cell"):
        return PROLIFERATIVE, "stem, progenitor, or mitotic germ-cell stage"
    if label in {"germline cell", "spermatocyte"}:
        return MEIOTIC_OR_MIXED, "mixed or meiotic germ-cell annotation"
    if label in {"spermatid"}:
        return POST_MITOTIC, "post-meiotic male germ cell"
    if _has(label, r"neuron|kenyon cell|photoreceptor|retinula|muscle|cardiomyocyte|glial|glia"):
        return POST_MITOTIC, "adult neural, glial, or muscle cell"
    if _has(label, r"fat body|fat mass|oenocyte|salivary gland|cone cell|pigment cell"):
        return POST_MITOTIC, "differentiated adult fly cell"
    if _has(label, r"enterocyte|enteroendocrine|copper cell|midgut large flat cell"):
        return POST_MITOTIC, "differentiated adult intestinal cell"
    if label in {"adult ventral nervous system", "adult peripheral nervous system"}:
        return POST_MITOTIC, "adult nervous-system compartment"
    if _has(label, r"hemocyte"):
        return PROLIFERATIVE, "adult fly immune cell with proliferative capacity"
    return UNCERTAIN, "mixed, generic, or insufficiently resolved annotation"


def classify_worm(label):
    if label == "Embroynic_cells":
        return PROLIFERATIVE, "dividing embryonic cells"
    if label == "Germline":
        return MEIOTIC_OR_MIXED, "mixed mitotic, meiotic, and developing germline"
    if label == "Sperms":
        return POST_MITOTIC, "terminally differentiated sperm"
    return POST_MITOTIC, "adult C. elegans somatic cell"


def classify_zebrafish(label):
    if _has(label, r"spermatogonia"):
        return PROLIFERATIVE, "mitotic spermatogonial stage"
    if _has(label, r"spermatocyte"):
        return MEIOTIC_OR_MIXED, "meiotic germ-cell stage"
    if _has(label, r"spermatid|spermatozoa"):
        return POST_MITOTIC, "post-meiotic male germ cell"
    if label in {"immune", "Sertoli & Leydig"}:
        return PROLIFERATIVE, "immune or somatic gonadal cell with proliferative capacity"
    return UNCERTAIN, "mixed, generic, or insufficiently resolved annotation"


def classify_cell_type(cell_type, species):
    label = str(cell_type).strip()
    lower = label.lower()
    if species in {"Mouse", "Human"}:
        group, basis = classify_mammal(lower)
    elif species == "Naked mole-rat":
        if lower == "basal":
            group, basis = PROLIFERATIVE, "epidermal basal progenitor compartment"
        elif lower == "corneum":
            group, basis = POST_MITOTIC, "terminally differentiated cornified compartment"
        else:
            group, basis = UNCERTAIN, "mixed or insufficiently resolved compartment"
    elif species == "Fruit fly":
        group, basis = classify_fly(lower)
    elif species == "C. elegans":
        group, basis = classify_worm(str(cell_type).strip())
    elif species == "Zebrafish":
        group, basis = classify_zebrafish(lower)
    else:
        raise ValueError(f"Unsupported species: {species}")
    return {"species": species, "cell_type": str(cell_type), "proliferation_class": group, "classification_basis": basis}
