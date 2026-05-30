import spaces
import gradio as gr
from gradio_molecule3d import Molecule3D
from gradio_cofoldinginput import CofoldingInput
import os
import re
import urllib.request
import yaml

# make sure to pip install gradio gradio_molecule3d gradio_cofoldinginput

CCD_URL = "https://huggingface.co/boltz-community/boltz-1/resolve/main/ccd.pkl"
MODEL_URL = "https://huggingface.co/boltz-community/boltz-1/resolve/main/boltz1.ckpt"

cache = "/home/user/.boltz"

os.makedirs(cache)

ccd = f"{cache}/ccd.pkl"
if not os.path.exists(ccd):
    print(
        f"Downloading the CCD dictionary to {ccd}. You may "
    )
    urllib.request.urlretrieve(CCD_URL, str(ccd))

# Download model
model =f"{cache}/boltz1.ckpt"
if not os.path.exists(model):
    print(
        f"Downloading the model weights to {model}"
    )
    urllib.request.urlretrieve(MODEL_URL, str(model))



@spaces.GPU(duration=120)
def predict(jobname, inputs, recycling_steps, sampling_steps, diffusion_samples):
    try:
        jobname = re.sub(r'[<>:"/\\|?*]', '_', jobname)
        if jobname == "":
            raise gr.Error("Job name empty or only invalid characters. Choose a plaintext name.")
        os.makedirs(jobname, exist_ok=True)
        """format Gradio Component:
        # {"chains": [
        #     {
        #         "class": "DNA",
        #         "sequence": "ATGCGT",
        #         "chain": "A"
        #     }
        # ], "covMods":[]
        # }
        """
        #sequences_for_msa = []
        output = {
        "sequences": []
        }
        representations = []
        for chain in inputs["chains"]:
            entity_type = chain["class"].lower()
            sequence_data = {
                entity_type: {
                    "id": chain["chain"],
                }
            }
            if entity_type in ["protein", "dna", "rna"]:
                sequence_data[entity_type]["sequence"] = chain["sequence"]
                if entity_type == "protein":
                    #sequences_for_msa.append(chain["sequence"])
                    if chain["msa"] == False:
                        sequence_data[entity_type]["msa"] = f"empty"
                representations.append({"model":0, "chain":chain["chain"], "style":"cartoon"})
            if entity_type == "ligand":
                if "sdf" in chain.keys():
                    if chain["sdf"]!="" and chain["name"]=="":
                        raise gr.Error("Sorry, no SDF support yet.")
                if "name" in chain.keys() and len(chain["name"])==3:
                     sequence_data[entity_type]["ccd"] = chain["name"]
                elif "smiles" in chain.keys():
                     sequence_data[entity_type]["smiles"] = chain["smiles"]
                else:
                    raise gr.Error("No ligand found, or not in the right format. CCD codes have 3 letters")
                    
                    
                representations.append({"model":0, "chain":chain["chain"], "style":"stick", "color":"greenCarbon"})
    
            if len(inputs["covMods"])>0:
                raise gr.Error("Sorry, covMods not supported yet. Coming soon. ")
            output["sequences"].append(sequence_data)
    
        # Convert the output to YAML
        yaml_file_path = f"{jobname}/{jobname}.yaml"
    
        # Write the YAML output to the file
        with open(yaml_file_path, "w") as file:
            yaml.dump(output, file, sort_keys=False, default_flow_style=False)
    
        os.system(f"cat {yaml_file_path}")
    
        os.system(f"boltz predict {jobname}/{jobname}.yaml --use_msa_server --out_dir {jobname} --recycling_steps {recycling_steps} --sampling_steps {sampling_steps} --diffusion_samples {diffusion_samples} --override --output_format pdb")
        print(os.listdir(jobname))
        print(os.listdir(f"{jobname}/boltz_results_{jobname}/predictions/{jobname}/"))
        return Molecule3D(f"{jobname}/boltz_results_{jobname}/predictions/{jobname}/{jobname}_model_0.pdb", label="Output", reps=representations)
    except Exception as e:
        raise gr.Error(f"failed with error:{e}")

with gr.Blocks() as blocks:
    gr.Markdown("# Boltz-1")
    gr.Markdown("""Open GUI for running [Boltz-1 model](https://github.com/jwohlwend/boltz/) <br>
    Key components:
    - MMSeqs2 Webserver [Mirdita et al.](https://www.nature.com/articles/s41592-022-01488-1)
    - Boltz-1 Model [Wohlwend et al.](https://github.com/jwohlwend/boltz/)
    - Gradio Custom Components [Molecule3D](https://huggingface.co/spaces/simonduerr/gradio_molecule3d)/[Cofolding Input](https://huggingface.co/spaces/simonduerr/gradio_cofoldinginput) by myself 
    - [3dmol.js Rego & Koes](https://academic.oup.com/bioinformatics/article/31/8/1322/213186)
    
    Note: This is an alpha: Some things like covalent modifications or using sdf files don't work yet. You can a Docker image of this on your local infrastructure easily using:
    `docker run -it -p 7860:7860 --platform=linux/amd64 --gpus all registry.hf.space/simonduerr-boltz-1:latest python app.py`
    """)
    with gr.Tab("Main"):
        jobname = gr.Textbox(label="Jobname")
        inp = CofoldingInput(label="Input")
        out = Molecule3D(label="Output")
    with gr.Tab("Settings"):
        recycling_steps =gr.Slider(value=3, minimum=0, label="Recycling steps")
        sampling_steps = gr.Slider(value=200, minimum=0, label="Sampling steps")
        diffusion_samples = gr.Slider(value=1, label="Diffusion samples")

    gr.Examples([
            ["TOP7",{"chains": [{"class": "protein", "msa":True,"sequence": "MGDIQVQVNIDDNGKNFDYTYTVTTESELQKVLNELMDYIKKQGAKRVRISITARTKKEAEKFAAILIKVFAELGYNDINVTFDGDTVTVEGQLEGGSLEHHHHHH","chain": "A"}], "covMods":[]}], 
            ["ApixacabanBinderSmiles", {"chains": [{"class": "protein", "msa":True,"sequence": "SVKSEYAEAAAVGQEAVAVFNTMKAAFQNGDKEAVAQYLARLASLYTRHEELLNRILEKARREGNKEAVTLMNEFTATFQTGKSIFNAMVAAFKNGDDDSFESYLQALEKVTAKGETLADQIAKAL","chain": "A"}, {"class":"ligand", "smiles":"COc1ccc(cc1)n2c3c(c(n2)C(=O)N)CCN(C3=O)c4ccc(cc4)N5CCCCC5=O", "sdf":"","name":"","chain": "B"}], "covMods":[]}], 
            ["ApixacabanBinderCCD", {"chains": [{"class": "protein","msa":True,"sequence": "SVKSEYAEAAAVGQEAVAVFNTMKAAFQNGDKEAVAQYLARLASLYTRHEELLNRILEKARREGNKEAVTLMNEFTATFQTGKSIFNAMVAAFKNGDDDSFESYLQALEKVTAKGETLADQIAKAL","chain": "A"}, {"class":"ligand", "name":"GG2", "sdf":"","chain": "B"}], "covMods":[]}] 
        ],
    inputs = [jobname, inp]     
    )

    btn = gr.Button("predict")

    btn.click(fn=predict, inputs=[jobname,inp, recycling_steps, sampling_steps, diffusion_samples], outputs=[out],  api_name="predict")

blocks.launch(ssr_mode=False)
