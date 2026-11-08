
import csv,re
ids=set("B-002,B-009,B-016,B-023,B-030,B-037,B-043,B-044,B-045,B-046,B-047,B-048,B-055,B-056,B-057,B-061,B-062,B-063,B-067,B-068,B-069,B-079,B-080,B-081,B-109,B-110,B-111,B-113,B-115,B-116,B-117,B-118,B-119,B-120,B-121,B-122,B-123,B-124,B-125,B-126,B-127,B-128,B-129,B-130,B-131,B-132,B-133,B-135,B-144,B-146,B-161,B-162,B-163,B-164,B-165,B-166,B-167,B-168,B-169,B-170,B-171,B-172,B-173,B-174,B-175,B-176,B-177,B-178,B-179,B-180,B-181,B-182,B-183,B-184,B-185,B-232,B-233,B-234,B-235,B-236,B-237,B-238,B-239,B-240,B-241,B-242,B-243,B-244,B-245,B-246,B-247,B-248,B-259,B-260,B-261,B-262,B-263,B-264,B-265,B-266,B-267,B-268,B-269,B-270,B-271,B-272,B-274,B-275,B-276,B-277,B-278,B-279,B-280,B-281,B-282,B-284,B-285,B-286,B-287,B-288,B-289,B-290,B-291,B-292,B-294,B-295,B-296,B-297,B-298,B-299,B-300,B-301,B-302,B-303,B-304,B-306,B-307,B-308,B-309,B-310,B-311,B-312,B-313,B-314,B-315,B-316,B-317,B-319,B-320,B-321,B-322,B-323,B-324,B-325,B-326,B-327,B-328,B-329,B-331,B-332,B-333,B-334,B-335,B-336,B-367,B-368,B-369,B-370,B-371,B-373,B-384,B-385,B-386,B-387,B-388,B-389,B-390,B-391,B-392,B-393,B-394,B-395,B-402,B-403,B-404,B-405,B-406,B-407,B-408,B-409,B-410,B-411,B-412,B-413,B-414,B-415,B-416,B-417,B-418,B-419".split(","))
rows=[]
for ln,line in enumerate(open("/home/user/2025.geo_spatialtrans/docs/BIO_REFERENCES.md",encoding="utf-8"),1):
    if not line.startswith("| B-"): continue
    cells=[x.strip() for x in line.rstrip("\n").strip().strip("|").split(" | ")]
    rows.append([ln]+cells)
allr=open("/tmp/r7_refs/bio_all.csv","w",newline="",encoding="utf-8"); w=csv.writer(allr)
w.writerow(["line","id","archetype","metric","value","interval","unit","source","session","notes","extra"])
n=0
for r in rows:
    r=r+[""]*(11-len(r)) if len(r)<11 else r[:10]+[" | ".join(r[10:])]
    w.writerow(r)
allr.close()
print(len(rows), sum(1 for r in rows if r[1] in ids), [i for i in ids if i not in {r[1] for r in rows}])
