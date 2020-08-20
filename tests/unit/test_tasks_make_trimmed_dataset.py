
import tempfile
import uuid

from pbcommand.testkit import PbIntegrationBase
from pbcommand.models.common import DataStore, DataStoreFile, FileTypes
from pbcore.io import ConsensusReadSet

import pbtestdata

class TestMakeTrimmedDataset(PbIntegrationBase):

    def test_integration(self):
        ccs_barcoded = pbtestdata.get_file("ccs-barcoded")
        datastore = tempfile.NamedTemporaryFile(suffix=".datastore.json").name
        lima_out = tempfile.NamedTemporaryFile(suffix=".consensusreadset.xml").name
        ccs_in = tempfile.NamedTemporaryFile(suffix=".consensusreadset.xml").name
        with ConsensusReadSet(ccs_barcoded) as ccs_tmp:
            ccs_tmp.name = "My Data (filtered)"
            ccs_tmp.tags = "ccs,filtered"
            ccs_tmp.write(ccs_in)
            ccs_tmp.name = "lima out"
            ccs_tmp.write(lima_out)
        ds = DataStore([
            DataStoreFile(uuid.uuid4(), "lima", FileTypes.DS_CCS.file_type_id, lima_out)
        ])
        ds.write_json(datastore)
        args = [
            "python3", "-m", "pbcoretools.tasks.make_trimmed_dataset",
            datastore, ccs_in
        ]
        self._check_call(args)
        with ConsensusReadSet("trimmed.consensusreadset.xml",
                              trustCounts=True) as ccs_out:
            assert ccs_out.numRecords > 0
            assert ccs_out.name == "My Data (trimmed)"
            assert ccs_out.tags == "ccs"
