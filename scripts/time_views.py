from argparse import ArgumentParser
from timeit import timeit
from urllib.parse import urljoin

import requests


class ViewTimer():

    object_types = ["cog", "pfam", "ko", "amr", "orthogroup"]

    @property
    def comparison_views(self):
        views = []
        views.extend([f"/{obj_type}_comparison" for obj_type in self.object_types])
        views.extend([f"/venn_{obj_type}/" for obj_type in self.object_types])
        views.extend([f"/pan_genome/{obj_type}/" for obj_type in self.object_types])
        views.extend([f"/plot_heatmap/{obj_type}/" for obj_type in self.object_types])
        views.extend([f"/entry_list_{obj_type}/" for obj_type in self.object_types])
        views.extend([f"/extract_{obj_type}/" for obj_type in self.object_types])
        return views

    def get_data(self, view, ntargets):
        if view == "/circos/":
            data = {"reference_taxid": "0",
                    "include_taxids": [str(i) for i in range(1, ntargets + 1)]}
            return data
        data = {"targets": [str(i) for i in range(ntargets)]}
        if view == "/amr_comparison":
            data["comp_type"] = "gene"
        return data

    def get_views(self):
        views = {"GET": ["orthogroup/group_0"],
                 "POST": self.comparison_views}
        views["POST"].append("/circos/")
        return views

    def time_views(self, base_url, number, ntargets):
        session = requests.session()

        # Get a csrftoken cookie and set the corresponding header
        session.get(urljoin(base_url, "venn_cog/"))
        session.headers.update({'X-CSRFToken': session.cookies['csrftoken']})
        results = []
        views = self.get_views()
        for view in views["GET"]:
            url = urljoin(base_url, view)
            res = timeit("session.get(url)",
                         globals={"session": session, "url": url},
                         number=number)
            results.append((view, res))
        for view in views["POST"]:
            url = urljoin(base_url, view)
            data = self.get_data(view, ntargets)
            res = timeit("session.post(url, data=data)",
                         globals={"session": session, "url": url, "data": data},
                         number=number)
            results.append((view, res))
        return results


if __name__ == "__main__":
    parser = ArgumentParser(prog="zDB")
    parser.add_argument("base_url")
    parser.add_argument("-f", "--filename", default="timings.csv")
    parser.add_argument("-n", "--number", default=5, type=int,
                        help="Number of times each view is called to compute "
                             "average rendering time")
    parser.add_argument("-t", "--targets", default=3, type=int,
                        help="Number of targets selected for comparison views")
    args = parser.parse_args()

    timer = ViewTimer()
    results = timer.time_views(args.base_url, args.number, args.targets)
    results = sorted(results, key=lambda x: x[1], reverse=True)

    with open(args.filename, "w") as fout:
        fout.write("\n".join(
            [",".join(map(str, [row[0], row[1] / args.number]))
             for row in results]))
