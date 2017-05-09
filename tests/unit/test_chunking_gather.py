
from collections import namedtuple
from zipfile import ZipFile
import tempfile
import tarfile
import shutil
import uuid
import json
import csv
import os
import unittest
import logging

from pbcommand.models.report import Report
from pbcommand.pb_io.report import load_report_from_json

import pbcoretools.chunking.gather as G

from base import get_temp_file, get_temp_dir

log = logging.getLogger(__name__)


class Record(object):
    def __init__(self, idx, alpha):
        self.idx = idx
        self.alpha = alpha

    def to_dict(self):
        return dict(id=self.idx, alpha=self.alpha)


def _to_n_records(nrecords):
    for i in xrange(nrecords):
        r = Record(i, 90)
        yield r


def _write_records_to_csv(records, output_csv):
    fields = records[0].to_dict().keys()
    with open(output_csv, 'w') as w:
        writer = csv.DictWriter(w, fieldnames=fields)
        writer.writeheader()
        writer.writerows([r.to_dict() for r in records])


class TestCsvGather(unittest.TestCase):

    def test_smoke(self):
        t = get_temp_file(suffix="-records-1.csv")
        _write_records_to_csv(list(_to_n_records(100)), t)

        t2 = get_temp_file(suffix="-records-2.csv")
        _write_records_to_csv(list(_to_n_records(57)), t2)

        tg = get_temp_file(suffix="records-gather.csv")
        G.gather_csv([t, t2], tg)

        nrecords = 0
        with open(tg, 'r') as r:
            reader = csv.DictReader(r)
            log.debug(reader.fieldnames)
            for _ in reader:
                nrecords += 1

        self.assertEqual(nrecords, 157)


def _write_stats_to_json(stats, output_json):
    with open(output_json, 'w') as w:
        w.write(Report.from_simple_dict(
            report_id="pbcommand_test",
            raw_d=stats,
            namespace="pb").to_json())


class TestJsonGather(unittest.TestCase):

    def test_smoke(self):
        t = get_temp_file(suffix="-stats-1.json")
        _write_stats_to_json({'n_reads':549,'n_zmws':100}, t)
        t2 = get_temp_file(suffix="-stats-2.json")
        _write_stats_to_json({'n_reads':733,'n_zmws':100}, t2)

        tg = get_temp_file(suffix="stats-gather.json")
        G.gather_report([t, t2], tg)

        r = load_report_from_json(tg)
        stats = { a.id:a.value for a in r.attributes }
        self.assertEqual(stats['pb_n_reads'], 549+733)
        self.assertEqual(stats['pb_n_zmws'], 200)



def _mkjson():
    uuid_ = str(uuid.uuid4())
    json_file = "{u}.json".format(u=uuid_)
    with open(json_file, "w") as json_out:
        json_out.write('{"uuid":"%s"}' % uuid_)
    return json_file


def create_tarball(file_name, n=2):
    with tarfile.open(file_name, mode="w:gz") as tgz_out:
        for i in range(n):
            json_file = _mkjson()
            tgz_out.add(json_file)
            os.remove(json_file)
    return file_name


class TestTgzGather(unittest.TestCase):

    def test_gather_tgz(self):
        tmp_dir = tempfile.mkdtemp()
        base_dir = os.getcwd()
        try:
            os.chdir(tmp_dir)
            tgz_files = [
                create_tarball(get_temp_file(suffix="-1.tgz")),
                create_tarball(get_temp_file(suffix="-2.tgz"))
            ]
            tgz_out = get_temp_file(suffix="-0.tgz")
            G.gather_tgz(tgz_files, tgz_out)
            uuids = set()
            with tarfile.open(tgz_out, mode="r:gz") as gathered:
                for member in gathered.getmembers():
                    d = json.loads(gathered.extractfile(member.name).read())
                    uuids.add(d["uuid"])
            self.assertEqual(len(uuids), 4)
        finally:
            os.chdir(base_dir)
            shutil.rmtree(tmp_dir)


def create_zip(file_name, n=2):
    with ZipFile(file_name, "w") as zip_out:
        for i in range(n):
            json_file = _mkjson()
            zip_out.write(json_file)
            os.remove(json_file)
    return file_name


class TestZipGather(unittest.TestCase):

    def test_gather_zip(self):
        tmp_dir = tempfile.mkdtemp()
        base_dir = os.getcwd()
        try:
            os.chdir(tmp_dir)
            zip_files = [
                create_zip(get_temp_file(suffix="-1.zip")),
                create_zip(get_temp_file(suffix="-2.zip"))
            ]
            zip_out = get_temp_file(suffix="-0.zip")
            G.gather_zip(zip_files, zip_out)
            uuids = set()
            with ZipFile(zip_out, "r") as gathered:
                for member in gathered.namelist():
                    d = json.loads(gathered.open(member).read())
                    uuids.add(d["uuid"])
            self.assertEqual(len(uuids), 4)
        finally:
            os.chdir(base_dir)
            shutil.rmtree(tmp_dir)
