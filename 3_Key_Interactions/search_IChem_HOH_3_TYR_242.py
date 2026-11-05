####
#REQUIRE PYTHON 3.12
####


import ichem_ifp

if __name__ == "__main__":

    ### IChem config ###

    config = ichem_ifp.IFPConfig()
    config.basic = True
    config.D_Ar = 5.5
    config.D_Hb = 3.6

    dico_interaction = {}

    protein = "protein_PGK2_cmp47.mol2" # To Change
    ligand = "ligand_PGK2_cmp47.mol2" # To Change

    ### IChem config ###

    ### Debut IFP ###

    lig_interactions = ichem_ifp.compute_ifp_interactions(
        protein, ligand, config
    )

    with open("interact_with_HOH_3_TYR_242.list", 'w') as file_output_hdl:
        for lig in lig_interactions:

            name = lig.ligand_name 

            interact_with_HOH_3 = False
            interact_with_TYR_242 = False

            for interaction in lig.interactions:
                type_interact = interaction.type_interaction
                residue = interaction.residue_identifier
                distance = interaction.distance
                atom_prot = interaction.atom_prot
                atom_lig = interaction.atom_lig
                angle = interaction.angle

                if type_interact.find("HBond_protein") != -1 and type_interact.find("Weak") == -1 and residue.find("HOH   3") != -1:
                    interact_with_HOH_3 = True
                if type_interact.find("HBond_protein") != -1 and type_interact.find("Weak") == -1 and residue.find("TYR 242") != -1:
                    interact_with_TYR_242 = True

            if interact_with_HOH_3 == True and interact_with_TYR_242 == True:
                file_output_hdl.write('{}\n'.format(name))


        ### Fin TC IFP ###
