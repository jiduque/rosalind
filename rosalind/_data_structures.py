from numbers import Number
from typing import Tuple, Optional, Generic, List


class Node(object):
    def __init__(self, key: Number, val: Generic) -> None:
        self.key = key
        self.val = val


class BinarySearchTree(object):
    def __init__(self, list_of_values: List[Tuple]) -> None:
        list_of_values.sort(key=lambda x: x[0])
        self.data = list(map(lambda x: Node(*x), list_of_values))
        self.n = len(self.data)

    def search(self, target: Number) -> Optional[Node]:
        left, right = 0, self.n - 1
        while left < right:
            mid = (left + right) // 2
            current_node = self.data[mid]
            if current_node.key == target:
                return current_node
            elif current_node.key < target:
                left = mid + 1
            else:
                right = mid - 1
        return None

    def search_closest(self, target: Number) -> Node:
        left, right = 0, self.n - 1
        if target >= self.data[right].key:
            return self.data[right]

        if target <= self.data[left].key:
            return self.data[left]

        while left < right:
            mid = (left + right) // 2
            current_node = self.data[mid]

            if current_node.key == target:
                return current_node

            if target < current_node.key:
                if mid > 0 and target > self.data[mid - 1].key:
                    return self._get_closest(current_node, self.data[mid - 1], target)
                right = mid

            else:
                if mid < self.n - 1 and target < self.data[mid + 1].key:
                    return self._get_closest(current_node, self.data[mid + 1], target)

                left = mid

        return current_node

    @staticmethod
    def _get_closest(x: Node, y: Node, target: Number):
        if abs(x.key - target) <= abs(y.key - target):
            return x
        return y


class UnionFind(object):
    def __init__(self, n: int) -> None:
        self.items = list(range(n))
        self.sizes = [1 for _ in range(n)]
        self.n_nodes = n
        self.n_connected_sets = n

    def find(self, p: int) -> int:
        root = self.items[p]
        while root != self.items[root]:
            self.items[root] = self.items[self.items[root]]
            root = self.items[root]
        self.items[p] = root
        return root

    def connected(self, p: int, q: int) -> bool:
        return p == q or self.find(p) == self.find(q)

    def union(self, p: int, q: int) -> None:
        if self.connected(p, q):
            return None

        p_root = self.find(p)
        q_root = self.find(q)

        smaller, larger = sorted([p_root, q_root], key=lambda z: self.sizes[z])
        self.items[smaller] = larger
        self.sizes[larger] += self.sizes[smaller]

        self.n_connected_sets -= 1

