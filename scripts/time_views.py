from argparse import ArgumentParser
from timeit import timeit
from urllib.parse import urljoin

import requests


class ViewTimer():

    object_types = ["cog", "pfam", "ko", "amr", "orthogroup"]

    @property
    def views(self):
        views = []
        views.extend([f"/{obj_type}_comparison" for obj_type in self.object_types])
        views.extend([f"/venn_{obj_type}/" for obj_type in self.object_types])
        views.extend([f"/pan_genome/{obj_type}/" for obj_type in self.object_types])
        views.extend([f"/plot_heatmap/{obj_type}/" for obj_type in self.object_types])
        views.extend([f"/entry_list_{obj_type}/" for obj_type in self.object_types])
        views.extend([f"/extract_{obj_type}/" for obj_type in self.object_types])
        return views

    def get_data(self, view, ntargets):
        data = {"targets": [str(i) for i in range(ntargets)]}
        if view == "/amr_comparison":
            data["comp_type"] = "gene"
        return data

    def time_views(self, base_url, number, ntargets):
        session = requests.session()

        # Get a csrftoken cookie and set the corresponding header
        session.get(urljoin(base_url, "venn_cog/"))
        session.headers.update({'X-CSRFToken': session.cookies['csrftoken']})
        results = []
        for view in self.views:
            url = urljoin(base_url, view)
            print(url)
            data = self.get_data(view, ntargets)
            res = timeit("session.post(url, data=data)",
                         globals={"session": session, "url": url, "data": data},
                         number=number)
            results.append((view, data, res))
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
    results = sorted(results, key=lambda x: x[2], reverse=True)

    with open(args.filename, "w") as fout:
        fout.write("\n".join(
            [",".join(map(str, [row[0], len(row[1]["targets"]), row[2] / args.number]))
             for row in results]))
