from qxalign import Qxalign
from memory_profiler import profile


@profile
def test_qxalign_simple():
    q = Qxalign()
    q.prepare("AAAACGT", "TGCA", "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    assert 60 == q.align()
    del q


@profile
def test_qxalign_long():
    q = Qxalign()

    q.prepare("AAAACGT", "TGCA", b"!!!!")
    assert 60 == q.align()
    q.trace()
    assert "3I 1=" == q.show_trace()

    q.prepare_query(query_seq="CAAC")
    assert 40 == q.align(semi=True)
    q.trace()
    assert "1X 3=" == q.show_trace()

    q.prepare("", "", "")
    try:
        q.align()
    except IndexError:
        pass
    else:
        assert False

    del q


if __name__ == "__main__":
    test_qxalign_simple()
