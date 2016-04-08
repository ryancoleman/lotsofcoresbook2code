def agts(queue):
    al = queue.add('surface.agts.py')
    queue.add('work_function.py', ncpus=1, deps=[al])

if __name__ == '__main__':
    execfile('Al100.py', {'k': 6, 'N': 5})
