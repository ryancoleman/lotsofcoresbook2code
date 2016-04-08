function open_row(row, id, tid) {
    request = new XMLHttpRequest();
    request.open('GET', '/open_row/' + id + '?x=' + tid, true);
    request.onload = function() {
        data = request.responseText;
        var table = document.getElementById('rows');
        r = row.rowIndex;
        if (data) {
            var newrow = table.insertRow(r + 1);
            var cell = newrow.insertCell(0);
            cell.colSpan = 100;
            cell.innerHTML = data;
        }
        else {
            table.deleteRow(r + 1);
        }
    }
    request.send();
}
