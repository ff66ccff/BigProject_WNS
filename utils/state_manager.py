import json
import os

class StateManager:
    def __init__(self, state_file="workflow_state.json"):
        self.state_file = state_file
        self.data = self._load_state()

    def _load_state(self):
        if os.path.exists(self.state_file):
            try:
                with open(self.state_file, 'r') as f:
                    return json.load(f)
            except json.JSONDecodeError:
                print(f"[Warning] State file {self.state_file} corrupted. Resetting.")
        # Default state
        return {
            "wrapper_cycle": 0,
            "shaker_stage": None,
            "current_receptor": "receptor.pdbqt"
        }

    def get(self, key, default=None):
        return self.data.get(key, default)

    def update(self, key, value):
        self.data[key] = value
        with open(self.state_file, 'w') as f:
            json.dump(self.data, f, indent=4)
        print(f"✅ [Checkpoint] State saved: {key} = {value}")