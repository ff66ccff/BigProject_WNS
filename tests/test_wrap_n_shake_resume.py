import unittest
import sys
import os
import json
from pathlib import Path
from unittest.mock import MagicMock, patch, mock_open

# Add project root and scripts folder to path
PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(PROJECT_ROOT))
sys.path.append(str(PROJECT_ROOT / "scripts"))

# Import the module to be tested
import wrap_n_shake_docking
from wrap_n_shake_docking import CheckpointState

class TestWrapNShakeResume(unittest.TestCase):
    def setUp(self):
        self.checkpoint_path = Path("fake_checkpoint.json")
        self.checkpoint = CheckpointState(self.checkpoint_path)

    def test_initial_state(self):
        """Test default state when initialized."""
        self.assertEqual(self.checkpoint.state["completed_seeds"], [])
        self.assertEqual(self.checkpoint.state["docked_ligand_files"], [])
        self.assertFalse(self.checkpoint.is_autogrid_complete())
        self.assertEqual(self.checkpoint.get_successful_docks(), 0)

    @patch("pathlib.Path.exists")
    @patch("pathlib.Path.open", new_callable=mock_open)
    def test_save_and_load(self, mock_file, mock_exists):
        """Test saving state to file and loading it back."""
        # Setup state
        self.checkpoint.mark_autogrid_complete()
        self.checkpoint.mark_seed_completed(101, "lig101.pdbqt", "rec.pdbqt")
        
        # Test Save
        self.checkpoint.save()
        mock_file.assert_called_with("w", encoding="utf-8")
        
        # Verify written JSON content
        handle = mock_file()
        args, _ = handle.write.call_args
        written_json = args[0]
        # We can't easily parse partial writes, but we can assume json.dump does its job if called.
        # Instead, let's mock the read back.
        
        # Mock Load
        mock_exists.return_value = True
        saved_state = {
            "completed_seeds": [101],
            "docked_ligand_files": ["lig101.pdbqt"],
            "current_receptor": "rec.pdbqt",
            "autogrid_complete": True,
            "successful_docks": 1
        }
        
        # Update mock open to read the saved state
        mock_file.side_effect = None # Reset side effect if any
        mock_file.return_value.__enter__.return_value.read.return_value = json.dumps(saved_state)
        
        # Create new instance to load
        new_checkpoint = CheckpointState(self.checkpoint_path)
        with patch("json.load", return_value=saved_state):
             loaded = new_checkpoint.load()
        
        self.assertTrue(loaded)
        self.assertTrue(new_checkpoint.is_autogrid_complete())
        self.assertTrue(new_checkpoint.is_seed_completed(101))
        self.assertEqual(new_checkpoint.get_current_receptor(), "rec.pdbqt")

    @patch("pathlib.Path.exists")
    def test_missing_checkpoint(self, mock_exists):
        """Test loading when checkpoint file is missing."""
        mock_exists.return_value = False
        loaded = self.checkpoint.load()
        self.assertFalse(loaded)
        # Should remain in default state
        self.assertEqual(self.checkpoint.state["completed_seeds"], [])

    @patch("pathlib.Path.exists")
    @patch("pathlib.Path.open", new_callable=mock_open, read_data="invalid json")
    def test_corrupt_checkpoint(self, mock_file, mock_exists):
        """Test loading a corrupt checkpoint file."""
        mock_exists.return_value = True
        # json.load will fail with the mock data, effectively simulating corruption if we don't mock json.load
        # but mock_open read_data is a string, json.load needs a file-like object.
        # Let's mock json.load to raise JSONDecodeError directly to be precise.
        
        with patch("json.load", side_effect=json.JSONDecodeError("msg", "doc", 0)):
            loaded = self.checkpoint.load()
            
        self.assertFalse(loaded)
        # Should reset/remain in default state
        self.assertEqual(self.checkpoint.state["completed_seeds"], [])

    def test_progression_logic(self):
        """Test the logic for marking seeds as complete."""
        # 1. Autogrid
        self.assertFalse(self.checkpoint.is_autogrid_complete())
        self.checkpoint.mark_autogrid_complete()
        self.assertTrue(self.checkpoint.is_autogrid_complete())
        
        # 2. First seed successful
        self.checkpoint.mark_seed_completed(101, "lig101.pdbqt", "rec_1.pdbqt")
        self.assertTrue(self.checkpoint.is_seed_completed(101))
        self.assertFalse(self.checkpoint.is_seed_completed(202))
        self.assertEqual(self.checkpoint.get_docked_ligands(), ["lig101.pdbqt"])
        self.assertEqual(self.checkpoint.get_current_receptor(), "rec_1.pdbqt")
        
        # 3. Second seed failed (no ligand file)
        self.checkpoint.mark_seed_completed(202, None, "rec_1.pdbqt")
        self.assertTrue(self.checkpoint.is_seed_completed(202))
        # Should not add to docked ligands list
        self.assertEqual(self.checkpoint.get_docked_ligands(), ["lig101.pdbqt"])
        # Successful docks count should not increase
        self.assertEqual(self.checkpoint.get_successful_docks(), 1)


if __name__ == "__main__":
    unittest.main()
