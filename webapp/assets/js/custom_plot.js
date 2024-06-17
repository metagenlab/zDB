// Note: dt is a DataTable with selected entries. The first column
// is expected to contain the entry identifier which can be passed
// to the custom_plots form.


function makeCustomPlot(dt) {
  let arr = [];
  let rows = dt.rows({selected: true})
  if (rows[0].length === 0) { rows = dt.rows() }
  rows.every(function(rowIdx, tableLoop, rowLoop) {
    // Extract data from columns 3 and 17 to put in JSON
    var data = this.data()
    arr.push(data[0]);
  });
  let entries = arr.join(',');
  let csrftoken = getCookie('csrftoken');
  $.redirect('/custom_plots/',
             {'entries': entries, 'csrfmiddlewaretoken': csrftoken},
             'POST',
             '_blank');
}
