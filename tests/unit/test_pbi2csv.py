
from pbcommand.testkit import PbIntegrationBase

from utils import skip_if_no_internal_data

class TestPbi2Csv(PbIntegrationBase):

    ALIGNMENTS = "/pbi/dept/secondary/siv/testdata/PacBioTestData/pbtestdata/data/AlignmentSet/m54075_160703_015941.alignmentset.xml"

    @skip_if_no_internal_data
    def test_pbi2csv(self):
        args = ["python", "-m", "pbcoretools.pbi2csv", self.ALIGNMENTS, "tst1.csv"]
        self._check_call(args)
        with open("tst1.csv") as csv_out:
            lines = csv_out.read().splitlines()
            self.assertEqual(lines[0:2], [
                "idx,file,hole,qstart,qend,qual,offset,flag,ref,tstart,tend,astart,aend,rc,matches,mismatches,mq,inserts,dels",
                "1,/pbi/dept/secondary/siv/testdata/PacBioTestData/pbtestdata/data/AlignmentSet/m54075_160703_015941.subreads.bam,34472412,12227,14732,0.8,55312384,3,lambda_NEB3011,13197,15634,12229,14732,TRUE,2212,122,254,169,103"
            ])

    @skip_if_no_internal_data
    def test_pbi2csv_load_snr(self):
        args = ["python", "-m", "pbcoretools.pbi2csv", self.ALIGNMENTS, "tst2.csv", "--load-snr"]
        self._check_call(args)
        with open("tst2.csv") as csv_out:
            lines = csv_out.read().splitlines()
            self.assertEqual(lines[0:2], [
                "idx,file,hole,qstart,qend,qual,offset,flag,ref,tstart,tend,astart,aend,rc,matches,mismatches,mq,inserts,dels,snrA,snrC,snrG,snrT",
                "1,/pbi/dept/secondary/siv/testdata/PacBioTestData/pbtestdata/data/AlignmentSet/m54075_160703_015941.subreads.bam,34472412,12227,14732,0.8,55312384,3,lambda_NEB3011,13197,15634,12229,14732,TRUE,2212,122,254,169,103,4.39776992798,8.43141269684,4.81357240677,7.40103912354"
            ])
