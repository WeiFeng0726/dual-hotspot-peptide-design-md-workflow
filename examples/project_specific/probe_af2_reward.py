import os
import traceback

import jax

from proteinfoundation.rewards.alphafold2_reward import AF2RewardModel
from proteinfoundation.utils.pdb_utils import get_chain_ids_from_pdb


PDB_PATH = (
    "/mnt/d/Project/cw/ly/Proteina-Complexa-runs/"
    "ly_b404_b405_pep8_12_probe/inference/"
    "search_binder_local_pipeline_LY_B404_B405_PEP8_12_ly_b404_b405_pep8_12_probe/"
    "job_0_n_167_id_0_single_orig0/job_0_n_167_id_0_single_orig0.pdb"
)
AF2_DIR = "/home/fatcat/Proteina-Complexa/community_models/ckpts/AF2"


def main() -> int:
    print("jax_version", jax.__version__)
    print("devices", jax.devices())
    print("default_backend", jax.default_backend())
    print(
        "env",
        {
            k: v
            for k, v in os.environ.items()
            if k.startswith(("JAX", "XLA", "CUDA", "NVIDIA"))
        },
    )

    target_chain, binder_chain = get_chain_ids_from_pdb(PDB_PATH)
    print("pdb_path", PDB_PATH)
    print("target_chain", target_chain)
    print("binder_chain", binder_chain)

    model = AF2RewardModel(
        protocol="binder",
        af_params_dir=AF2_DIR,
        reward_weights={
            "con": 0.0,
            "i_pae": -1.0,
            "plddt": 0.0,
            "dgram_cce": 0.0,
            "min_ipae": 0.0,
            "min_ipsae": 0.0,
            "avg_ipsae": 0.0,
            "max_ipsae": 0.0,
            "min_ipsae_10": 0.0,
            "max_ipsae_10": 0.0,
            "avg_ipsae_10": 0.0,
        },
        use_multimer=True,
        num_recycles=1,
        use_initial_guess=True,
        use_initial_atom_pos=False,
        seed=0,
        device_id=0,
    )
    try:
        result = model.score(
            PDB_PATH,
            requires_grad=False,
            binder_chain=binder_chain,
            target_chain=target_chain,
        )
        print("score_keys", sorted(result.keys()))
        print("total_reward", result["total_reward"])
        print("reward", result["reward"])
        return 0
    except Exception as exc:
        print("EXC_TYPE", type(exc).__name__)
        print("EXC", exc)
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
